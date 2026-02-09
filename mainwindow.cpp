#include "mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    // 1. Setup UI
    setupGlobalUI();

    // 2. Setup Audio Engine
    m_player = new QMediaPlayer(this);
    m_audioOutput = new QAudioOutput(this);
    m_player->setAudioOutput(m_audioOutput);
    m_audioOutput->setVolume(1.0);

    // 3. Setup Timer
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &MainWindow::updateTracker);
}

MainWindow::~MainWindow()
{
}

// =========================================================
// UI SETUP
// =========================================================

void MainWindow::setupGlobalUI()
{
    resize(1200, 800);
    setWindowTitle("Test of Bassline  - Adjusting loop start and stop crashes");

    // Global Menu
    QMenu *fileMenu = menuBar()->addMenu("File");
    QAction *actOpen = new QAction("Open WAV...", this);
    actOpen->setShortcut(QKeySequence::Open);
    connect(actOpen, &QAction::triggered, this, &MainWindow::openFile);
    fileMenu->addAction(actOpen);

    QAction *actExit = new QAction("Exit", this);
    connect(actExit, &QAction::triggered, this, &QWidget::close);
    fileMenu->addAction(actExit);

    // Main Tab Holder
    m_tabs = new QTabWidget(this);
    setCentralWidget(m_tabs);

    // Init All Tabs
    initAnalyserTab();      // 1. Note Grid
    initSpectrogramTab();   // 2. Heatmap
    initEnergyTab();        // 3. Rhythm

}

// --- TAB 1: STANDARD analyser ---
void MainWindow::initAnalyserTab()
{
    QWidget *tab = new QWidget();
    QVBoxLayout *mainLayout = new QVBoxLayout(tab);

    // Controls
    QHBoxLayout *ctrlLayout = new QHBoxLayout();
    m_btnPlay = new QPushButton("▶ Play");
    m_btnStop = new QPushButton("⏹ Stop");
    m_thresholdLabel = new QLabel("Signal Threshold: 1.00");
    m_thresholdSlider = new QSlider(Qt::Horizontal);
    m_thresholdSlider->setRange(1, 300);
    m_thresholdSlider->setValue(100);
    connect(m_thresholdSlider, &QSlider::valueChanged, this, &MainWindow::onThresholdChanged);

    ctrlLayout->addWidget(m_btnPlay);
    ctrlLayout->addWidget(m_btnStop);
    ctrlLayout->addWidget(m_thresholdLabel);
    ctrlLayout->addWidget(m_thresholdSlider);
    ctrlLayout->addStretch();
    mainLayout->addLayout(ctrlLayout);

    // Loop Controls
    QHBoxLayout *loopLayout = new QHBoxLayout();
    m_lblLoopStart = new QLabel("Loop Start (0%)");
    m_loopStartSlider = new QSlider(Qt::Horizontal);
    m_loopStartSlider->setRange(0, 1000);
    m_loopStartSlider->setValue(0);
    m_lblLoopEnd = new QLabel("Loop End (100%)");
    m_loopEndSlider = new QSlider(Qt::Horizontal);
    m_loopEndSlider->setRange(0, 1000);
    m_loopEndSlider->setValue(1000);
    connect(m_loopStartSlider, &QSlider::valueChanged, this, &MainWindow::onLoopSliderChanged);
    connect(m_loopEndSlider, &QSlider::valueChanged, this, &MainWindow::onLoopSliderChanged);

    loopLayout->addWidget(m_lblLoopStart);
    loopLayout->addWidget(m_loopStartSlider);
    loopLayout->addWidget(m_lblLoopEnd);
    loopLayout->addWidget(m_loopEndSlider);
    mainLayout->addLayout(loopLayout);

    // Splitter
    QSplitter *splitter = new QSplitter(Qt::Horizontal);
    m_table = new QTableWidget();
    m_table->setColumnCount(3);
    m_table->setHorizontalHeaderLabels({"Step", "Freq(Hz)", "Note"});
    splitter->addWidget(m_table);

    m_plot = new QCustomPlot();
    m_plot->addGraph();
    m_plot->graph(0)->setPen(QPen(Qt::black));
    m_plot->xAxis->setLabel("Time (s)");
    m_plot->yAxis->setLabel("Amplitude");
    m_plot->setInteraction(QCP::iRangeDrag, true);
    m_plot->setInteraction(QCP::iRangeZoom, true);

    m_lineStart = new QCPItemLine(m_plot); m_lineStart->setPen(QPen(Qt::green, 2));
    m_lineEnd = new QCPItemLine(m_plot);   m_lineEnd->setPen(QPen(Qt::red, 2));

    splitter->addWidget(m_plot);
    splitter->setStretchFactor(1, 4);
    mainLayout->addWidget(splitter);

    connect(m_btnPlay, &QPushButton::clicked, this, &MainWindow::startPlayback);
    connect(m_btnStop, &QPushButton::clicked, this, &MainWindow::stopPlayback);

    m_tabs->addTab(tab, "1. Waveform & Notes");
}

// --- TAB 2: SPECTROGRAM ---
void MainWindow::initSpectrogramTab()
{
    QWidget *tab = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(tab);

    QLabel *info = new QLabel("Spectrogram: X = Time, Y = Frequency, Color = Loudness.");
    layout->addWidget(info);

    m_spectrogramPlot = new QCustomPlot();
    m_spectrogram = new QCPColorMap(m_spectrogramPlot->xAxis, m_spectrogramPlot->yAxis);
    m_colorScale = new QCPColorScale(m_spectrogramPlot);
    m_spectrogramPlot->plotLayout()->addElement(0, 1, m_colorScale);
    m_spectrogram->setColorScale(m_colorScale);
    QCPColorGradient gradient;
    gradient.loadPreset(QCPColorGradient::gpPolar);
    m_spectrogram->setGradient(gradient);
    m_spectrogram->setInterpolate(false);
    m_spectrogramPlot->xAxis->setLabel("Time (s)");
    m_spectrogramPlot->yAxis->setLabel("Frequency (Hz)");
    m_spectrogramPlot->setInteraction(QCP::iRangeDrag, true);
    m_spectrogramPlot->setInteraction(QCP::iRangeZoom, true);

    layout->addWidget(m_spectrogramPlot);
    m_tabs->addTab(tab, "2. Spectrogram");
}

// --- TAB 3: ENERGY ---
void MainWindow::initEnergyTab()
{
    QWidget *tab = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(tab);

    QHBoxLayout *top = new QHBoxLayout();
    m_bpmLabel = new QLabel("Estimated BPM: --");
    QPushButton *btnCalc = new QPushButton("Recalculate Rhythm");
    connect(btnCalc, &QPushButton::clicked, this, &MainWindow::calculateBPM);
    top->addWidget(m_bpmLabel);
    top->addWidget(btnCalc);
    top->addStretch();
    layout->addLayout(top);

    m_energyPlot = new QCustomPlot();
    m_energyPlot->addGraph();
    m_energyPlot->graph(0)->setPen(QPen(Qt::blue, 2));
    m_energyPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20)));
    m_energyPlot->xAxis->setLabel("Time (s)");
    m_energyPlot->yAxis->setLabel("Energy (RMS)");
    m_energyPlot->setInteraction(QCP::iRangeDrag, true);
    m_energyPlot->setInteraction(QCP::iRangeZoom, true);

    layout->addWidget(m_energyPlot);
    m_tabs->addTab(tab, "3. Rhythm Energy");
}



// =========================================================
// LOGIC
// =========================================================

void MainWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Open WAV", "", "WAV Files (*.wav)");
    if (fileName.isEmpty()) return;

    m_currentFileName = fileName;

    bool ok;
    int steps = QInputDialog::getInt(this, "Steps", "Analysis Steps:", 16, 8, 64, 4, &ok);
    if (!ok) return;
    m_numSteps = steps;

    if (loadWavFile(fileName)) {
        // 1. Basic Plot
        int nSamples = m_audioData.size();
        QVector<double> x(nSamples), y(nSamples);
        for (int i=0; i<nSamples; ++i) {
            x[i] = (double)i / m_sampleRate;
            y[i] = m_audioData[i];
        }
        m_plot->graph(0)->setData(x, y);
        m_plot->rescaleAxes();

        // Reset State
        m_loopStartSlider->setValue(0);
        m_loopEndSlider->setValue(1000);
        m_loopStartSample = 0;
        m_loopEndSample = nSamples;
        double totalTime = (double)nSamples / m_sampleRate;

        m_lineStart->start->setCoords(0, -1); m_lineStart->end->setCoords(0, 1);
        m_lineEnd->start->setCoords(totalTime, -1); m_lineEnd->end->setCoords(totalTime, 1);

        updateVisualSteps();

        m_plot->replot();

        m_player->setSource(QUrl::fromLocalFile(fileName));

        // RUN analysers
        analyseSteps();         // Tab 1
        generateSpectrogram();  // Tab 2
        analyseEnergy();        // Tab 3

        statusBar()->showMessage("Analysis Complete: " + fileName);
    }
}




// --- STANDARD FUNCTIONS ---

void MainWindow::onThresholdChanged(int value)
{
    m_magnitudeThreshold = value / 100.0;
    m_thresholdLabel->setText(QString("Signal Threshold: %1").arg(m_magnitudeThreshold, 0, 'f', 2));
    analyseSteps();
}

void MainWindow::onLoopSliderChanged()
{
    if (m_audioData.empty()) return;
    int startVal = m_loopStartSlider->value();
    int endVal = m_loopEndSlider->value();
    int minGap = 10;

    if (sender() == m_loopStartSlider) {
        if (startVal > endVal - minGap) {
            m_loopEndSlider->blockSignals(true);
            m_loopEndSlider->setValue(qMin(1000, startVal + minGap));
            m_loopEndSlider->blockSignals(false);
            endVal = m_loopEndSlider->value();
        }
    } else {
        if (endVal < startVal + minGap) {
            m_loopStartSlider->blockSignals(true);
            m_loopStartSlider->setValue(qMax(0, endVal - minGap));
            m_loopStartSlider->blockSignals(false);
            startVal = m_loopStartSlider->value();
        }
    }

    int totalSamples = m_audioData.size();
    m_loopStartSample = (startVal / 1000.0) * totalSamples;
    m_loopEndSample   = (endVal / 1000.0) * totalSamples;

    m_lblLoopStart->setText(QString("Loop Start (%1%)").arg(startVal / 10.0));
    m_lblLoopEnd->setText(QString("Loop End (%1%)").arg(endVal / 10.0));

    double duration = (double)totalSamples / m_sampleRate;
    double startTime = (startVal / 1000.0) * duration;
    double endTime = (endVal / 1000.0) * duration;

    m_lineStart->start->setCoords(startTime, -1); m_lineStart->end->setCoords(startTime, 1);
    m_lineEnd->start->setCoords(endTime, -1); m_lineEnd->end->setCoords(endTime, 1);

    updateVisualSteps();
    analyseSteps();
    m_plot->replot();
}

void MainWindow::updateVisualSteps()
{
    for (auto line : m_stepLines) { if (line) { m_plot->removeItem(line); delete line; } }
    m_stepLines.clear();

    if (m_audioData.empty()) return;
    double startSec = (double)m_loopStartSample / m_sampleRate;
    double endSec = (double)m_loopEndSample / m_sampleRate;
    double loopDuration = endSec - startSec;
    if (loopDuration <= 0.001) return;
    double stepTime = loopDuration / m_numSteps;

    for (int i = 1; i < m_numSteps; ++i) {
        double pos = startSec + (i * stepTime);
        QCPItemLine *line = new QCPItemLine(m_plot);
        line->start->setCoords(pos, -1);
        line->end->setCoords(pos, 1);
        line->setPen(QPen(Qt::red, 1, Qt::DashLine));
        m_stepLines.push_back(line);
    }
}

void MainWindow::startPlayback()
{
    if (m_currentFileName.isEmpty()) { openFile(); return; }
    double startSec = (double)m_loopStartSample / m_sampleRate;
    m_player->setPosition(startSec * 1000);
    m_player->play();
    m_timer->start(30);
}

void MainWindow::stopPlayback()
{
    m_player->stop();
    m_timer->stop();
    m_table->clearSelection();
}

void MainWindow::updateTracker()
{
    if (m_player->playbackState() != QMediaPlayer::PlayingState) {
        if (m_player->playbackState() == QMediaPlayer::StoppedState) stopPlayback();
        return;
    }
    qint64 posMs = m_player->position();
    double posSec = posMs / 1000.0;
    double loopStartSec = (double)m_loopStartSample / m_sampleRate;
    double loopEndSec   = (double)m_loopEndSample / m_sampleRate;

    if (posSec >= loopEndSec || posSec < loopStartSec) {
        m_player->setPosition(loopStartSec * 1000);
        return;
    }

    double loopDuration = loopEndSec - loopStartSec;
    if (loopDuration > 0.001) {
        double currentPosInLoop = posSec - loopStartSec;
        double progress = currentPosInLoop / loopDuration;
        int currentStep = floor(progress * m_numSteps);
        if (currentStep >= 0 && currentStep < m_numSteps) m_table->selectRow(currentStep);
    }
}

void MainWindow::analyseSteps()
{
    if (m_audioData.empty()) return;
    int activeRegionSize = m_loopEndSample - m_loopStartSample;
    if (activeRegionSize < 64) return;

    m_stepSize = (double)activeRegionSize / m_numSteps;
    m_table->setRowCount(m_numSteps);

    const int N_FFT = 8192;
    kiss_fft_cfg cfg = kiss_fft_alloc(N_FFT, 0, NULL, NULL);
    kiss_fft_cpx in[N_FFT], out[N_FFT];

    for (int i=0; i<m_numSteps; ++i) {
        int startIdx = m_loopStartSample + floor(i * m_stepSize);
        for (int k=0; k<N_FFT; ++k) {
            if (k < m_stepSize && (startIdx + k) < (int)m_audioData.size()) {
                float win = 1.0f;
                if (m_stepSize > 2.0) win = 0.5 * (1 - cos(2 * M_PI * k / (m_stepSize - 1)));
                in[k].r = m_audioData[startIdx + k] * win;
                in[k].i = 0;
            } else {
                in[k].r = 0; in[k].i = 0;
            }
        }
        kiss_fft(cfg, in, out);
        float maxMag = 0; int maxIdx = 0;
        int startBin = 30 * N_FFT / m_sampleRate;
        int endBin = 500 * N_FFT / m_sampleRate;

        for (int k=startBin; k<endBin; ++k) {
            float mag = sqrt(out[k].r*out[k].r + out[k].i*out[k].i);
            if (mag > maxMag) { maxMag = mag; maxIdx = k; }
        }
        float domFreq = maxIdx * m_sampleRate / (float)N_FFT;
        m_table->setItem(i, 0, new QTableWidgetItem(QString::number(i+1)));
        if (maxMag > m_magnitudeThreshold) {
            m_table->setItem(i, 1, new QTableWidgetItem(QString::number(domFreq, 'f', 1)));
            m_table->setItem(i, 2, new QTableWidgetItem(freqToNote(domFreq)));
        } else {
            m_table->setItem(i, 1, new QTableWidgetItem("-"));
            m_table->setItem(i, 2, new QTableWidgetItem("Quiet"));
        }
    }
    free(cfg);
}

void MainWindow::generateSpectrogram()
{
    if (m_audioData.empty()) return;
    int nSamples = m_audioData.size();
    int fftSize = 1024; int overlap = 512;
    int timeSteps = (nSamples - fftSize) / (fftSize - overlap);
    if (timeSteps <= 0) return;
    int freqBins = fftSize / 2;

    m_spectrogram->data()->setSize(timeSteps, freqBins);
    m_spectrogram->data()->setRange(QCPRange(0, (double)nSamples/m_sampleRate), QCPRange(0, m_sampleRate/2));

    kiss_fft_cfg cfg = kiss_fft_alloc(fftSize, 0, NULL, NULL);
    kiss_fft_cpx in[fftSize], out[fftSize];

    for (int t = 0; t < timeSteps; ++t) {
        int startSample = t * (fftSize - overlap);
        for (int i = 0; i < fftSize; ++i) {
            if (startSample + i < nSamples) {
                float win = 0.5 * (1 - cos(2 * M_PI * i / (fftSize - 1)));
                in[i].r = m_audioData[startSample + i] * win; in[i].i = 0;
            } else { in[i].r = 0; in[i].i = 0; }
        }
        kiss_fft(cfg, in, out);
        for (int f = 0; f < freqBins; ++f) {
            double mag = sqrt(out[f].r * out[f].r + out[f].i * out[f].i);
            double db = 20 * log10(mag + 1e-6);
            m_spectrogram->data()->setCell(t, f, db);
        }
    }
    free(cfg);
    m_spectrogramPlot->rescaleAxes();
    m_spectrogramPlot->replot();
}

void MainWindow::analyseEnergy()
{
    if (m_audioData.empty()) return;
    int windowSize = m_sampleRate / 50;
    int nSamples = m_audioData.size();
    int steps = nSamples / windowSize;

    QVector<double> x(steps), y(steps);
    for (int i = 0; i < steps; ++i) {
        double sumSq = 0;
        int start = i * windowSize;
        for (int j = 0; j < windowSize; ++j) {
            if (start + j < nSamples) {
                float val = m_audioData[start + j];
                sumSq += val * val;
            }
        }
        x[i] = (double)start / m_sampleRate;
        y[i] = sqrt(sumSq / windowSize);
    }
    m_energyPlot->graph(0)->setData(x, y);
    m_energyPlot->rescaleAxes();
    m_energyPlot->replot();
    calculateBPM();
}

void MainWindow::calculateBPM()
{
    auto data = m_energyPlot->graph(0)->data();
    if (data->isEmpty()) return;
    double maxVal = 0;
    for (auto it = data->begin(); it != data->end(); ++it) if (it->value > maxVal) maxVal = it->value;
    double threshold = maxVal * 0.6;
    int beatCount = 0; bool inBeat = false;
    double firstBeatTime = -1; double lastBeatTime = -1;

    for (auto it = data->begin(); it != data->end(); ++it) {
        if (!inBeat && it->value > threshold) {
            inBeat = true; beatCount++;
            if (firstBeatTime < 0) firstBeatTime = it->key;
            lastBeatTime = it->key;
        } else if (inBeat && it->value < threshold * 0.8) { inBeat = false; }
    }

    if (beatCount < 2) { m_bpmLabel->setText("BPM: ??"); return; }
    double duration = lastBeatTime - firstBeatTime;
    if (duration > 0) {
        double bpm = ((beatCount - 1) / duration) * 60.0;
        m_bpmLabel->setText(QString("Estimated BPM: %1").arg(QString::number(bpm, 'f', 1)));
    }
}

QString MainWindow::freqToNote(float freq)
{
    if (freq <= 0) return "-";
    const QString names[] = {"C","C#","D","D#","E","F","F#","G","G#","A","A#","B"};
    float C0 = 16.35;
    int h = qRound(12 * log2(freq / C0));
    int octave = h / 12;
    int note = h % 12;
    if (note < 0) note = 0; if (note > 11) note = 11;
    return names[note] + QString::number(octave);
}

bool MainWindow::loadWavFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly)) return false;
    struct RiffHeader { char chunkID[4]; uint32_t chunkSize; char format[4]; } riff;
    if (file.read((char*)&riff, sizeof(RiffHeader)) != sizeof(RiffHeader)) return false;

    bool foundFmt = false, foundData = false;
    uint16_t numChannels = 0;

    while (!file.atEnd()) {
        char chunkID[4]; uint32_t chunkSize;
        if (file.read(chunkID, 4) != 4) break;
        if (file.read((char*)&chunkSize, 4) != 4) break;

        if (qstrncmp(chunkID, "fmt ", 4) == 0) {
            struct FmtChunk {
                uint16_t audioFormat; uint16_t numChannels; uint32_t sampleRate;
                uint32_t byteRate; uint16_t blockAlign; uint16_t bitsPerSample;
            } fmt;
            if (file.read((char*)&fmt, 16) != 16) break;
            m_sampleRate = fmt.sampleRate; m_bitsPerSample = fmt.bitsPerSample;
            numChannels = fmt.numChannels; foundFmt = true;
            if (chunkSize > 16) file.seek(file.pos() + chunkSize - 16);
        } else if (qstrncmp(chunkID, "data", 4) == 0) {
            if (!foundFmt) break;
            QByteArray raw = file.read(chunkSize);
            int bytesPerSample = m_bitsPerSample / 8;
            if (bytesPerSample == 0) break;
            int frames = (raw.size() / bytesPerSample) / numChannels;
            m_audioData.resize(frames);
            const char* ptr = raw.constData();
            for (int i = 0; i < frames; ++i) {
                float sum = 0.0f;
                for (int c = 0; c < numChannels; ++c) {
                    float val = 0.0f;
                    if (m_bitsPerSample == 16) { val = (*reinterpret_cast<const int16_t*>(ptr)) / 32768.0f; ptr += 2; }
                    else if (m_bitsPerSample == 24) { int32_t v = 0; v |= (uint8_t)ptr[0]; v |= (uint8_t)ptr[1]<<8; v |= (int8_t)ptr[2]<<16; val = v/8388608.0f; ptr+=3; }
                    else if (m_bitsPerSample == 32) { val = *reinterpret_cast<const float*>(ptr); ptr += 4; }
                    else if (m_bitsPerSample == 8) { val = (*reinterpret_cast<const quint8*>(ptr) - 128) / 128.0f; ptr += 1; }
                    sum += val;
                }
                m_audioData[i] = sum / numChannels;
            }
            foundData = true; break;
        } else { file.seek(file.pos() + chunkSize); }
    }
    return foundData;
}
