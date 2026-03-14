#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTabWidget>
#include <QTableWidget>
#include <QPushButton>
#include <QMenu>
#include <QAction>
#include <QMenuBar>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSplitter>
#include <QLabel>
#include <QSlider>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <QMediaPlayer>
#include <QAudioOutput>
#include <QTimer>
#include <QtMath>
#include <vector>

// Libs
#include "libs/qcustomplot.h"
extern "C" {
#include "libs/kiss_fft.h"
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    // Shared Slots
    void openFile();
    void startPlayback();
    void stopPlayback();
    void updateTracker();
    void runFullAnalysis();

    void onBassMaxFreqChanged(int value);

    void onThresholdChanged(int value);
    void onLoopSliderChanged();
    void calculateBPM();
    void analyseChords();

    void exportBasslineLMMS();
    void exportChordLMMS();

private:
    QTabWidget *m_tabs;
    QTimer *m_analysisTimer;

    // --- TAB 1: ANALYSER ---
    QCustomPlot *m_plot;
    QTableWidget *m_table;
    QPushButton *m_btnPlay;
    QPushButton *m_btnStop;
    QSlider *m_thresholdSlider;
    QLabel *m_thresholdLabel;

    QSlider *m_loopStartSlider;
    QSlider *m_loopEndSlider;
    QLabel *m_lblLoopStart;
    QLabel *m_lblLoopEnd;
    QCPItemLine *m_lineStart;
    QCPItemLine *m_lineEnd;
    std::vector<QCPItemLine*> m_stepLines;
    QSlider *m_bassMaxFreqSlider;
    QLabel *m_bassMaxFreqLabel;
    int m_bassMaxFreq = 500;

    // --- TAB 2: SPECTROGRAM ---
    QCustomPlot *m_spectrogramPlot;
    QCPColorMap *m_spectrogram;
    QCPColorScale *m_colorScale;

    // --- TAB 3: ENERGY / BPM ---
    QCustomPlot *m_energyPlot;
    QLabel *m_bpmLabel;

    // --- TAB 4: CHORD ---
    QString identifyChord(const std::vector<int>& midiNotes);
    int freqToMidi(float freq);
    QSpinBox *m_chordMinFreq;
    QSpinBox *m_chordMaxFreq;
    QTableWidget *m_chordTable;
    QCustomPlot *m_chordSpectrumPlot;
    QSpinBox *m_notchFreq;
    QSpinBox *m_notchWidth;

    // --- DATA STATE ---
    int m_loopStartSample = 0;
    int m_loopEndSample = 0;
    std::vector<float> m_audioData;
    int m_sampleRate = 44100;
    int m_numSteps = 16;
    double m_stepSize = 0;
    int m_bitsPerSample = 16;
    double m_magnitudeThreshold = 1.0;
    QString m_currentFileName;

    // Audio Engine
    QMediaPlayer *m_player;
    QAudioOutput *m_audioOutput;
    QTimer *m_timer;

    // --- INIT FUNCTIONS ---
    void setupGlobalUI();
    void initAnalyserTab();
    void initSpectrogramTab();
    void initEnergyTab();
    void initChordTab();


    // --- HELPERS ---
    void analyseSteps();
    void updateVisualSteps();
    void generateSpectrogram();
    void analyseEnergy();

    QString freqToNote(float freq);
    bool loadWavFile(const QString &fileName);

    std::vector<int> m_basslineMidi;
    std::vector<std::vector<int>> m_chordMidi;
};

#endif // MAINWINDOW_H
