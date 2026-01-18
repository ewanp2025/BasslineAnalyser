#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMediaPlayer>
#include <QAudioOutput>
#include <QTimer>
#include <vector>
#include <QUrl>

// Libraries
#include "libs/qcustomplot.h"
extern "C" {
#include "libs/kiss_fft.h"
}

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_actionOpen_triggered(); // File -> Open
    void startPlayback();
    void stopPlayback();
    void updateTracker();

private:
    Ui::MainWindow *ui;

    // Audio Data
    std::vector<float> m_audioData;
    int m_sampleRate = 44100;
    int m_numSteps = 16;
    int m_stepSize = 0;

    // File Info
    QString m_currentFileName;
    int m_bitsPerSample = 16;

    // The New "Bulletproof" Player
    QMediaPlayer *m_player;
    QAudioOutput *m_audioOutput;
    QTimer *m_timer;

    // Helper functions
    void analyzeSteps();
    QString freqToNote(float freq);
    bool loadWavFile(const QString &fileName);
};
#endif // MAINWINDOW_H
