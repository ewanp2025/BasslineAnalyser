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

    // UI Slots
    void onThresholdChanged(int value);
    void onLoopSliderChanged();
    void calculateBPM(); // Tab 3

private:
    // --- UI ELEMENTS ---
    QTabWidget *m_tabs;

    // --- TAB 1: ANALYSER ---
    QCustomPlot *m_plot;
    QTableWidget *m_table;
    QPushButton *m_btnPlay;
    QPushButton *m_btnStop;
    QSlider *m_thresholdSlider;
    QLabel *m_thresholdLabel;

    // Loop Controls
    QSlider *m_loopStartSlider;
    QSlider *m_loopEndSlider;
    QLabel *m_lblLoopStart;
    QLabel *m_lblLoopEnd;
    QCPItemLine *m_lineStart;
    QCPItemLine *m_lineEnd;
    std::vector<QCPItemLine*> m_stepLines;

    // --- TAB 2: SPECTROGRAM ---
    QCustomPlot *m_spectrogramPlot;
    QCPColorMap *m_spectrogram;
    QCPColorScale *m_colorScale;

    // --- TAB 3: ENERGY / BPM ---
    QCustomPlot *m_energyPlot;
    QLabel *m_bpmLabel;

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


    // --- HELPERS ---
    void analyseSteps();
    void updateVisualSteps();
    void generateSpectrogram();
    void analyseEnergy();

    QString freqToNote(float freq);
    bool loadWavFile(const QString &fileName);
};

#endif // MAINWINDOW_H
