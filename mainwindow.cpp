#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <QtMath>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // 1. Setup Plot
    ui->plot->addGraph();
    ui->plot->graph(0)->setPen(QPen(Qt::black));
    ui->plot->xAxis->setLabel("Time (s)");
    ui->plot->yAxis->setLabel("Amplitude");
    ui->plot->setInteraction(QCP::iRangeDrag, true);
    ui->plot->setInteraction(QCP::iRangeZoom, true);

    // 2. Setup Table
    QStringList headers;
    headers << "Step" << "Freq(Hz)" << "Note";
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setHorizontalHeaderLabels(headers);
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    // 3. Setup The "Bulletproof" Player
    m_player = new QMediaPlayer(this);
    m_audioOutput = new QAudioOutput(this);
    m_player->setAudioOutput(m_audioOutput);
    m_audioOutput->setVolume(1.0); // Full volume

    // 4. Timer for Tracker
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &MainWindow::updateTracker);

    // 5. Connect Buttons
    connect(ui->btnPlay, &QPushButton::clicked, this, &MainWindow::startPlayback);
    connect(ui->btnStop, &QPushButton::clicked, this, &MainWindow::stopPlayback);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionOpen_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open WAV", "", "WAV Files (*.wav)");
    if (fileName.isEmpty()) return;

    m_currentFileName = fileName;

    bool ok;
    int steps = QInputDialog::getInt(this, "Steps", "Number of steps (16, 32, 64):", 16, 16, 64, 16, &ok);
    if (!ok) return;
    m_numSteps = steps;

    if (loadWavFile(fileName)) {
        // Plot Waveform
        QVector<double> x(m_audioData.size()), y(m_audioData.size());
        for (size_t i=0; i<m_audioData.size(); ++i) {
            x[i] = (double)i / m_sampleRate;
            y[i] = m_audioData[i];
        }
        ui->plot->graph(0)->setData(x, y);
        ui->plot->rescaleAxes();

        // Draw Grid Lines
        ui->plot->clearItems();
        double totalTime = (double)m_audioData.size() / m_sampleRate;
        double stepTime = totalTime / m_numSteps;
        for(int i=1; i<m_numSteps; ++i) {
            QCPItemLine *line = new QCPItemLine(ui->plot);
            line->start->setCoords(i*stepTime, -1);
            line->end->setCoords(i*stepTime, 1);
            line->setPen(QPen(Qt::red, 1, Qt::DashLine));
        }
        ui->plot->replot();

        // Run Analysis
        analyzeSteps();

        // Setup Player Source
        m_player->setSource(QUrl::fromLocalFile(fileName));
    }
}

bool MainWindow::loadWavFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly)) return false;

    struct RiffHeader { char chunkID[4]; uint32_t chunkSize; char format[4]; } riff;
    file.read((char*)&riff, sizeof(RiffHeader));

    bool foundData = false;
    uint16_t numChannels = 1;

    while (!file.atEnd()) {
        char chunkID[4]; uint32_t chunkSize;
        if (file.read(chunkID, 4) != 4) break;
        if (file.read((char*)&chunkSize, 4) != 4) break;

        if (qstrncmp(chunkID, "fmt ", 4) == 0) {
            struct FmtChunk {
                uint16_t audioFormat; uint16_t numChannels; uint32_t sampleRate;
                uint32_t byteRate; uint16_t blockAlign; uint16_t bitsPerSample;
            } fmt;
            file.read((char*)&fmt, sizeof(FmtChunk));
            m_sampleRate = fmt.sampleRate;
            m_bitsPerSample = fmt.bitsPerSample;
            numChannels = fmt.numChannels;
            if (chunkSize > sizeof(FmtChunk)) file.seek(file.pos() + chunkSize - sizeof(FmtChunk));
        }
        else if (qstrncmp(chunkID, "data", 4) == 0) {
            QByteArray raw = file.read(chunkSize);

            int bytesPerSample = m_bitsPerSample / 8;
            int totalSamples = raw.size() / bytesPerSample;
            int frames = totalSamples / numChannels;

            m_audioData.resize(frames);
            const char* ptr = raw.constData();

            for (int i=0; i<frames; ++i) {
                float sum = 0.0f;
                for (int c=0; c<numChannels; ++c) {
                    float sampleVal = 0.0f;
                    if (m_bitsPerSample == 16) {
                        sampleVal = (*reinterpret_cast<const int16_t*>(ptr)) / 32768.0f;
                        ptr += 2;
                    } else if (m_bitsPerSample == 24) {
                        int32_t val = 0;
                        val |= (uint8_t)(ptr[0]) << 0;
                        val |= (uint8_t)(ptr[1]) << 8;
                        val |= (int8_t) (ptr[2]) << 16;
                        sampleVal = val / 8388608.0f;
                        ptr += 3;
                    } else if (m_bitsPerSample == 32) {
                        sampleVal = *reinterpret_cast<const float*>(ptr);
                        ptr += 4;
                    }
                    sum += sampleVal;
                }
                m_audioData[i] = sum / numChannels;
            }
            foundData = true;
            break;
        }
        else {
            file.seek(file.pos() + chunkSize);
        }
    }
    return foundData;
}

void MainWindow::analyzeSteps()
{
    if (m_audioData.empty()) return;

    m_stepSize = m_audioData.size() / m_numSteps;
    ui->tableWidget->setRowCount(m_numSteps);

    const int N_FFT = 8192;
    kiss_fft_cfg cfg = kiss_fft_alloc(N_FFT, 0, NULL, NULL);
    kiss_fft_cpx in[N_FFT], out[N_FFT];

    qDebug() << "--- START ANALYSIS ---";

    for (int i=0; i<m_numSteps; ++i) {
        int startIdx = i * m_stepSize;

        // 1. Fill Buffer
        for (int k=0; k<N_FFT; ++k) {
            if (k < m_stepSize && (startIdx + k) < m_audioData.size()) {
                float win = 0.5 * (1 - cos(2 * M_PI * k / (m_stepSize - 1)));
                in[k].r = m_audioData[startIdx + k] * win;
                in[k].i = 0;
            } else {
                in[k].r = 0; in[k].i = 0;
            }
        }

        kiss_fft(cfg, in, out);

        // 2. Find Peak
        float maxMag = 0;
        int maxIdx = 0;
        int startBin = 30 * N_FFT / m_sampleRate;
        int endBin = 500 * N_FFT / m_sampleRate;

        for (int k=startBin; k<endBin; ++k) {
            float mag = sqrt(out[k].r*out[k].r + out[k].i*out[k].i);
            if (mag > maxMag) {
                maxMag = mag;
                maxIdx = k;
            }
        }

        float domFreq = maxIdx * m_sampleRate / (float)N_FFT;

        // DEBUG: Print the first 5 steps to see what's happening
        if (i < 5) qDebug() << "Step" << i << "Mag:" << maxMag << "Freq:" << domFreq;

        // 3. Populate Table
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(QString::number(i+1)));

        //THRESHOLD
        if (maxMag > 1.0) {
            ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(domFreq, 'f', 1)));
            ui->tableWidget->setItem(i, 2, new QTableWidgetItem(freqToNote(domFreq)));
        } else {
            ui->tableWidget->setItem(i, 1, new QTableWidgetItem("-"));
            ui->tableWidget->setItem(i, 2, new QTableWidgetItem("Quiet"));
        }
    }
    free(cfg);
}

QString MainWindow::freqToNote(float freq)
{
    if (freq <= 0) return "-";
    const QString names[] = {"C","C#","D","D#","E","F","F#","G","G#","A","A#","B"};
    float C0 = 16.35;
    int h = qRound(12 * log2(freq / C0));
    int octave = h / 12;
    int note = h % 12;
    return names[note] + QString::number(octave);
}

void MainWindow::startPlayback()
{
    if (m_currentFileName.isEmpty()) return;

    // Simple, reliable playback
    m_player->stop();
    m_player->play();
    m_timer->start(50);
}

void MainWindow::stopPlayback()
{
    m_player->stop();
    m_timer->stop();
    ui->tableWidget->clearSelection();
}

void MainWindow::updateTracker()
{
    if (m_player->playbackState() != QMediaPlayer::PlayingState) {
        if (m_player->playbackState() == QMediaPlayer::StoppedState) stopPlayback();
        return;
    }

    qint64 pos = m_player->position();
    qint64 dur = m_player->duration();

    if (dur > 0) {
        double progress = (double)pos / dur;
        int currentStep = floor(progress * m_numSteps);

        if (currentStep < m_numSteps) {
            ui->tableWidget->selectRow(currentStep);
        }
    }
}
