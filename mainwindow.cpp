#include "mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{

    setupGlobalUI();

    m_player = new QMediaPlayer(this);
    m_audioOutput = new QAudioOutput(this);
    m_player->setAudioOutput(m_audioOutput);
    m_audioOutput->setVolume(1.0);

    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &MainWindow::updateTracker);
    m_analysisTimer = new QTimer(this);
    m_analysisTimer->setSingleShot(true); // Ensures it only fires once after a delay
    connect(m_analysisTimer, &QTimer::timeout, this, &MainWindow::runFullAnalysis);
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
    setWindowTitle("Bassline & Chord Analyser");

    QMenu *fileMenu = menuBar()->addMenu("File");
    QAction *actOpen = new QAction("Open WAV...", this);
    actOpen->setShortcut(QKeySequence::Open);
    connect(actOpen, &QAction::triggered, this, &MainWindow::openFile);
    fileMenu->addAction(actOpen);

    QAction *actExit = new QAction("Exit", this);
    connect(actExit, &QAction::triggered, this, &QWidget::close);
    fileMenu->addAction(actExit);


    QWidget *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    QVBoxLayout *masterLayout = new QVBoxLayout(centralWidget);

    QHBoxLayout *globalCtrlLayout = new QHBoxLayout();
    m_btnPlay = new QPushButton("▶ Play");
    m_btnStop = new QPushButton("⏹ Stop");
    connect(m_btnPlay, &QPushButton::clicked, this, &MainWindow::startPlayback);
    connect(m_btnStop, &QPushButton::clicked, this, &MainWindow::stopPlayback);

    globalCtrlLayout->addWidget(m_btnPlay);
    globalCtrlLayout->addWidget(m_btnStop);
    globalCtrlLayout->addStretch();
    masterLayout->addLayout(globalCtrlLayout);

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
    masterLayout->addLayout(loopLayout);

    m_plot = new QCustomPlot();
    m_plot->setMinimumHeight(200); // Restrict height so it doesn't take over the screen
    m_plot->setMaximumHeight(250);
    m_plot->addGraph();
    m_plot->graph(0)->setPen(QPen(Qt::black));
    m_plot->xAxis->setLabel("Time (s)");
    m_plot->yAxis->setLabel("Amplitude");
    m_plot->setInteraction(QCP::iRangeDrag, true);
    m_plot->setInteraction(QCP::iRangeZoom, true);

    m_lineStart = new QCPItemLine(m_plot); m_lineStart->setPen(QPen(Qt::green, 2));
    m_lineEnd = new QCPItemLine(m_plot);   m_lineEnd->setPen(QPen(Qt::red, 2));

    masterLayout->addWidget(m_plot);

    m_tabs = new QTabWidget(this);
    masterLayout->addWidget(m_tabs);

    initAnalyserTab();
    initSpectrogramTab();
    initEnergyTab();
    initChordTab();
}

// --- TAB 1: Baseline analyser ---
void MainWindow::initAnalyserTab()
{
    QWidget *tab = new QWidget();
    QVBoxLayout *mainLayout = new QVBoxLayout(tab);

    QHBoxLayout *ctrlLayout = new QHBoxLayout();
    m_thresholdLabel = new QLabel("Signal Threshold: 1.00");
    m_thresholdSlider = new QSlider(Qt::Horizontal);
    m_thresholdSlider->setRange(1, 300);
    m_thresholdSlider->setValue(100);
    connect(m_thresholdSlider, &QSlider::valueChanged, this, &MainWindow::onThresholdChanged);

    m_bassMaxFreqLabel = new QLabel("Max Bass Freq: 500 Hz");
    m_bassMaxFreqSlider = new QSlider(Qt::Horizontal);
    m_bassMaxFreqSlider->setRange(60, 2000);
    m_bassMaxFreqSlider->setValue(500);
    connect(m_bassMaxFreqSlider, &QSlider::valueChanged, this, &MainWindow::onBassMaxFreqChanged);

    QPushButton *btnExportBass = new QPushButton("Export LMMS");
    connect(btnExportBass, &QPushButton::clicked, this, &MainWindow::exportBasslineLMMS);

    ctrlLayout->addWidget(m_thresholdLabel);
    ctrlLayout->addWidget(m_thresholdSlider);
    ctrlLayout->addWidget(m_bassMaxFreqLabel);
    ctrlLayout->addWidget(m_bassMaxFreqSlider);
    ctrlLayout->addWidget(btnExportBass);
    ctrlLayout->addStretch();
    mainLayout->addLayout(ctrlLayout);

    m_table = new QTableWidget();
    m_table->setColumnCount(3);
    m_table->setHorizontalHeaderLabels({"Step", "Freq(Hz)", "Note"});
    mainLayout->addWidget(m_table);

    m_tabs->addTab(tab, "1. Bassline Notes");
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

// --- TAB 4: CHORD ---
void MainWindow::initChordTab()
{
    QWidget *tab = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(tab);

    // Controls & Notch Filter
    QHBoxLayout *ctrlLayout = new QHBoxLayout();

    ctrlLayout->addWidget(new QLabel("Min Freq:"));
    m_chordMinFreq = new QSpinBox(); m_chordMinFreq->setRange(20, 5000); m_chordMinFreq->setValue(200);
    ctrlLayout->addWidget(m_chordMinFreq);

    ctrlLayout->addWidget(new QLabel("Max Freq:"));
    m_chordMaxFreq = new QSpinBox(); m_chordMaxFreq->setRange(100, 10000); m_chordMaxFreq->setValue(2000);
    ctrlLayout->addWidget(m_chordMaxFreq);

    ctrlLayout->addWidget(new QLabel(" | Notch Freq (0=Off):"));
    m_notchFreq = new QSpinBox(); m_notchFreq->setRange(0, 5000); m_notchFreq->setValue(0);
    ctrlLayout->addWidget(m_notchFreq);

    ctrlLayout->addWidget(new QLabel("Notch Width:"));
    m_notchWidth = new QSpinBox(); m_notchWidth->setRange(10, 500); m_notchWidth->setValue(50);
    ctrlLayout->addWidget(m_notchWidth);

    QPushButton *btnAnalyse = new QPushButton("Analyse Chords");
    connect(btnAnalyse, &QPushButton::clicked, this, &MainWindow::analyseChords);

    QPushButton *btnExportChord = new QPushButton("Export LMMS");
    connect(btnExportChord, &QPushButton::clicked, this, &MainWindow::exportChordLMMS);

    ctrlLayout->addWidget(btnAnalyse);
    ctrlLayout->addWidget(btnExportChord);
    ctrlLayout->addStretch();
    layout->addLayout(ctrlLayout);

    QSplitter *splitter = new QSplitter(Qt::Vertical);

    m_chordSpectrumPlot = new QCustomPlot();
    m_chordSpectrumPlot->addGraph();
    m_chordSpectrumPlot->graph(0)->setPen(QPen(Qt::blue));
    m_chordSpectrumPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 50)));
    m_chordSpectrumPlot->xAxis->setLabel("Frequency (Hz)");
    m_chordSpectrumPlot->yAxis->setLabel("Average Magnitude");
    splitter->addWidget(m_chordSpectrumPlot);

    m_chordTable = new QTableWidget();
    m_chordTable->setColumnCount(4);
    m_chordTable->setHorizontalHeaderLabels({"Step", "Raw Notes", "Pitch Classes", "Detected Chord"});
    m_chordTable->horizontalHeader()->setStretchLastSection(true);
    splitter->addWidget(m_chordTable);

    splitter->setStretchFactor(0, 2); // Give the plot more space than the table
    layout->addWidget(splitter);

    m_tabs->addTab(tab, "4. Chord Analyser");
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
        analyseChords();

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
    m_plot->replot();

    updateVisualSteps();
    m_analysisTimer->start(150);
    m_plot->replot();
}

void MainWindow::updateVisualSteps()
{
    if (m_audioData.empty()) return;
    double startSec = (double)m_loopStartSample / m_sampleRate;
    double endSec = (double)m_loopEndSample / m_sampleRate;
    double loopDuration = endSec - startSec;
    if (loopDuration <= 0.001) return;
    double stepTime = loopDuration / m_numSteps;


    while (m_stepLines.size() < (size_t)(m_numSteps - 1)) {
        QCPItemLine *line = new QCPItemLine(m_plot);
        line->setPen(QPen(Qt::red, 1, Qt::DashLine));
        m_stepLines.push_back(line);
    }


    for (size_t i = 0; i < m_stepLines.size(); ++i) {
        if (i < (size_t)(m_numSteps - 1)) {
            // Move the line to the new position
            double pos = startSec + ((i + 1) * stepTime);
            m_stepLines[i]->start->setCoords(pos, -1);
            m_stepLines[i]->end->setCoords(pos, 1);
            m_stepLines[i]->setVisible(true); // Ensure it is shown
        } else {
            m_stepLines[i]->setVisible(false);
        }
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

    m_basslineMidi.assign(m_numSteps, -1); // Reset with -1 (no note)

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
        int endBin = m_bassMaxFreq * N_FFT / m_sampleRate;
        endBin = std::min(N_FFT / 2 - 1, endBin); // Safety clamp

        for (int k=startBin; k<endBin; ++k) {
            float mag = sqrt(out[k].r*out[k].r + out[k].i*out[k].i);
            if (mag > maxMag) { maxMag = mag; maxIdx = k; }
        }
        float domFreq = maxIdx * m_sampleRate / (float)N_FFT;
        m_table->setItem(i, 0, new QTableWidgetItem(QString::number(i+1)));

        if (maxMag > m_magnitudeThreshold) {
            m_table->setItem(i, 1, new QTableWidgetItem(QString::number(domFreq, 'f', 1)));
            m_table->setItem(i, 2, new QTableWidgetItem(freqToNote(domFreq)));
            m_basslineMidi[i] = freqToMidi(domFreq);
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
    const int fftSize = 1024; 
    const int overlap = 512;
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
    if (data->isEmpty() || data->size() < 2) return;

    int n = data->size();
    std::vector<double> energy(n);
    int i = 0;
    for (auto it = data->begin(); it != data->end(); ++it) {
        energy[i++] = it->value;
    }

    std::vector<double> autocorr(n, 0.0);
    for (int lag = 0; lag < n / 2; ++lag) { // Only need to check up to half the window
        double sum = 0;
        for (int t = 0; t < n - lag; ++t) {
            sum += energy[t] * energy[t + lag];
        }
        autocorr[lag] = sum;
    }


    double dt = (data->at(1)->key - data->at(0)->key); // time between energy steps
    int minLag = std::ceil((60.0 / 200.0) / dt); // Min lag for 200 BPM max
    int maxLag = std::min((int)(n / 2), (int)std::ceil((60.0 / 60.0) / dt)); // Max lag for 60 BPM min

    double maxCorr = 0;
    int bestLag = -1;

    for (int lag = minLag; lag < maxLag; ++lag) {
        // Simple peak picking: must be greater than neighbors
        if (autocorr[lag] > autocorr[lag - 1] && autocorr[lag] > autocorr[lag + 1]) {
            if (autocorr[lag] > maxCorr) {
                maxCorr = autocorr[lag];
                bestLag = lag;
            }
        }
    }


    if (bestLag > 0) {
        double beatDuration = bestLag * dt;
        double bpm = 60.0 / beatDuration;
        m_bpmLabel->setText(QString("Estimated BPM (Autocorr): %1").arg(bpm, 0, 'f', 1));
    } else {
        m_bpmLabel->setText("Estimated BPM: ??");
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


int MainWindow::freqToMidi(float freq) {
    if (freq <= 0) return 0;
    return qRound(12 * log2(freq / 440.0) + 69);
}

void MainWindow::analyseChords()
{
    if (m_audioData.empty()) return;
    int activeRegionSize = m_loopEndSample - m_loopStartSample;
    if (activeRegionSize < 64) return;

    m_stepSize = (double)activeRegionSize / m_numSteps;
    m_chordTable->setRowCount(m_numSteps);

    const int N_FFT = 8192;
    kiss_fft_cfg cfg = kiss_fft_alloc(N_FFT, 0, NULL, NULL);
    kiss_fft_cpx in[N_FFT], out[N_FFT];

    int minBin = m_chordMinFreq->value() * N_FFT / m_sampleRate;
    int maxBin = m_chordMaxFreq->value() * N_FFT / m_sampleRate;
    minBin = std::max(1, minBin);
    maxBin = std::min(N_FFT / 2 - 1, maxBin);


    int notchCenter = m_notchFreq->value();
    int notchW = m_notchWidth->value();
    int notchMinBin = std::max(0, (notchCenter - notchW/2) * N_FFT / m_sampleRate);
    int notchMaxBin = std::min(N_FFT/2 - 1, (notchCenter + notchW/2) * N_FFT / m_sampleRate);
    bool useNotch = (notchCenter > 0);


    std::vector<double> avgSpectrum(N_FFT / 2, 0.0);

    m_chordMidi.assign(m_numSteps, std::vector<int>()); // Reset

    for (int i = 0; i < m_numSteps; ++i) {
        int startIdx = m_loopStartSample + floor(i * m_stepSize);

        for (int k = 0; k < N_FFT; ++k) {
            if (k < m_stepSize && (startIdx + k) < (int)m_audioData.size()) {
                float win = 0.5 * (1 - cos(2 * M_PI * k / (m_stepSize - 1)));
                in[k].r = m_audioData[startIdx + k] * win;
                in[k].i = 0;
            } else {
                in[k].r = 0; in[k].i = 0;
            }
        }

        kiss_fft(cfg, in, out);


        for (int k = 0; k < N_FFT / 2; ++k) {
            float mag = sqrt(out[k].r*out[k].r + out[k].i*out[k].i);
            avgSpectrum[k] += mag / m_numSteps;
        }

        struct Peak { int bin; float mag; };
        std::vector<Peak> peaks;

        for (int k = minBin; k < maxBin; ++k) {

            if (useNotch && k >= notchMinBin && k <= notchMaxBin) continue;

            float mag = sqrt(out[k].r*out[k].r + out[k].i*out[k].i);
            float prevMag = sqrt(out[k-1].r*out[k-1].r + out[k-1].i*out[k-1].i);
            float nextMag = sqrt(out[k+1].r*out[k+1].r + out[k+1].i*out[k+1].i);

            if (mag > prevMag && mag > nextMag && mag > m_magnitudeThreshold * 0.5) {
                peaks.push_back({k, mag});
            }
        }

        std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
            return a.mag > b.mag;
        });

        std::vector<int> midiNotes;
        QString noteStr = "";
        for (size_t p = 0; p < std::min((size_t)4, peaks.size()); ++p) {
            float freq = peaks[p].bin * m_sampleRate / (float)N_FFT;
            int midi = freqToMidi(freq);
            if (std::find(midiNotes.begin(), midiNotes.end(), midi) == midiNotes.end()) {
                midiNotes.push_back(midi);
                noteStr += freqToNote(freq) + " ";
            }
        }

        m_chordMidi[i] = midiNotes;

        QString chordName = identifyChord(midiNotes);
        QString pcStr = "";
        for(int n : midiNotes) pcStr += QString::number(n % 12) + " ";

        m_chordTable->setItem(i, 0, new QTableWidgetItem(QString::number(i + 1)));
        m_chordTable->setItem(i, 1, new QTableWidgetItem(noteStr.trimmed()));
        m_chordTable->setItem(i, 2, new QTableWidgetItem(pcStr.trimmed()));
        m_chordTable->setItem(i, 3, new QTableWidgetItem(chordName));
    }
    free(cfg);


    QVector<double> specX(N_FFT / 2), specY(N_FFT / 2);
    for (int k = 0; k < N_FFT / 2; ++k) {
        specX[k] = k * m_sampleRate / (double)N_FFT;
        specY[k] = avgSpectrum[k];
    }
    m_chordSpectrumPlot->graph(0)->setData(specX, specY);
    m_chordSpectrumPlot->rescaleAxes();
    m_chordSpectrumPlot->xAxis->setRange(m_chordMinFreq->value() * 0.8, m_chordMaxFreq->value() * 1.2);
    m_chordSpectrumPlot->replot();
}


QString MainWindow::identifyChord(const std::vector<int>& midiNotes)
{
    if (midiNotes.size() < 2) return "Single Note / Unknown";


    std::vector<int> sortedNotes = midiNotes;
    std::sort(sortedNotes.begin(), sortedNotes.end());

    int root = sortedNotes[0];
    std::vector<int> intervals;
    for (int note : sortedNotes) {
        int interval = (note - root) % 12;
        if (std::find(intervals.begin(), intervals.end(), interval) == intervals.end()) {
            intervals.push_back(interval);
        }
    }
    std::sort(intervals.begin(), intervals.end());


    bool hasMinorThird = std::find(intervals.begin(), intervals.end(), 3) != intervals.end();
    bool hasMajorThird = std::find(intervals.begin(), intervals.end(), 4) != intervals.end();
    bool hasFifth = std::find(intervals.begin(), intervals.end(), 7) != intervals.end();

    QString rootName = freqToNote(440.0 * pow(2.0, (root - 69) / 12.0));
    rootName.chop(1);

    if (hasMajorThird && hasFifth) return rootName + " Major";
    if (hasMinorThird && hasFifth) return rootName + " Minor";
    if (hasMajorThird && !hasFifth) return rootName + " Major (implied)";
    if (hasMinorThird && !hasFifth) return rootName + " Minor (implied)";
    if (std::find(intervals.begin(), intervals.end(), 7) != intervals.end() && !hasMajorThird && !hasMinorThird)
        return rootName + " Power Chord (5th)";

    return rootName + " (Complex/Unknown)";
}

void MainWindow::onBassMaxFreqChanged(int value)
{
    m_bassMaxFreq = value;
    m_bassMaxFreqLabel->setText(QString("Max Bass Freq: %1 Hz").arg(value));
    analyseSteps();
}

void MainWindow::runFullAnalysis()
{
    analyseSteps();
    analyseChords();
}


void MainWindow::exportBasslineLMMS()
{
    if (m_basslineMidi.empty()) {
        QMessageBox::warning(this, "Export Error", "No bassline data to export. Open a file first.");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "Save Bassline Pattern", "bassline.xpt", "LMMS Pattern (*.xpt)");
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;

    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<!DOCTYPE lmms-project>\n";
    out << "<lmms-project creator=\"BasslineAnalyser\" version=\"20\" type=\"pattern\">\n";
    out << "  <head/>\n";
    out << "  <pattern muted=\"0\" name=\"Detected Bass\" type=\"1\" steps=\"" << m_numSteps << "\" pos=\"0\">\n";

    for (int i = 0; i < m_numSteps; ++i) {
        if (m_basslineMidi[i] != -1) {
            // Write out the note. 'pos' advances by 48 ticks per step. 'len' is 48 (1 step long).
            out << "    <note key=\"" << m_basslineMidi[i] << "\" pan=\"0\" vol=\"100\" len=\"48\" pos=\"" << (i * 48) << "\"/>\n";
        }
    }

    out << "  </pattern>\n";
    out << "</lmms-project>\n";

    file.close();
    statusBar()->showMessage("Exported Bassline to " + fileName, 3000);
}

void MainWindow::exportChordLMMS()
{
    if (m_chordMidi.empty()) {
        QMessageBox::warning(this, "Export Error", "No chord data to export. Open a file first.");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "Save Chord Pattern", "chords.xpt", "LMMS Pattern (*.xpt)");
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;

    QTextStream out(&file);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<!DOCTYPE lmms-project>\n";
    out << "<lmms-project creator=\"BasslineAnalyser\" version=\"20\" type=\"pattern\">\n";
    out << "  <head/>\n";
    out << "  <pattern muted=\"0\" name=\"Detected Chords\" type=\"1\" steps=\"" << m_numSteps << "\" pos=\"0\">\n";

    for (int i = 0; i < m_numSteps; ++i) {
        for (int note : m_chordMidi[i]) {
            // Write out every note detected in the chord for this step
            out << "    <note key=\"" << note << "\" pan=\"0\" vol=\"100\" len=\"48\" pos=\"" << (i * 48) << "\"/>\n";
        }
    }

    out << "  </pattern>\n";
    out << "</lmms-project>\n";

    file.close();
    statusBar()->showMessage("Exported Chords to " + fileName, 3000);
}
