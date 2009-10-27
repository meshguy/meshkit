/********************************************************************************
** Form generated from reading ui file 'itaps.ui'
**
** Created: Mon May 4 17:24:30 2009
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_ITAPS_H
#define UI_ITAPS_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QAction *actionExit;
    QWidget *centralwidget;
    QWidget *widget_DisplayArea;
    QFrame *frame;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBox_displayGeometricNodes;
    QCheckBox *checkBox_DisplayContinuousGeomCurves;
    QCheckBox *checkBox_DisplayDiscretizedGeomCurves;
    QLabel *label_3;
    QSpinBox *spinBox_CurrentEdgeID;
    QLabel *label_8;
    QLineEdit *lineEdit_NumOfEdgeSegments;
    QLabel *label_7;
    QWidget *gridLayoutWidget_2;
    QGridLayout *gridLayout_2;
    QLineEdit *lineEdit_MinEdgeLength;
    QLabel *label_5;
    QLineEdit *lineEdit_MaxEdgeLength;
    QLabel *label_4;
    QLabel *label_6;
    QLineEdit *lineEdit_MeanEdgeLength;
    QFrame *line_2;
    QLabel *label;
    QRadioButton *radioButton_DisplayUVSurface;
    QRadioButton *radioButton_DisplayXYZSurface;
    QLabel *label_2;
    QSpinBox *spinBox_CurrentFaceID;
    QPushButton *pushButton_GenerateSurface;
    QPushButton *pushButton_AcceptSurface;
    QPushButton *pushButton_RejectSurface;
    QLabel *label_9;
    QCheckBox *checkBox_displayCurrentEdge;
    QCheckBox *checkBox_DisplayCurrentFace;
    QFrame *line;
    QMenuBar *menubar;
    QMenu *menuFile;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1151, 838);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        widget_DisplayArea = new QWidget(centralwidget);
        widget_DisplayArea->setObjectName(QString::fromUtf8("widget_DisplayArea"));
        widget_DisplayArea->setGeometry(QRect(10, 10, 851, 771));
        frame = new QFrame(centralwidget);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setGeometry(QRect(870, 10, 271, 621));
        frame->setFrameShape(QFrame::StyledPanel);
        frame->setFrameShadow(QFrame::Raised);
        verticalLayoutWidget = new QWidget(frame);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 30, 231, 92));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        checkBox_displayGeometricNodes = new QCheckBox(verticalLayoutWidget);
        checkBox_displayGeometricNodes->setObjectName(QString::fromUtf8("checkBox_displayGeometricNodes"));

        verticalLayout->addWidget(checkBox_displayGeometricNodes);

        checkBox_DisplayContinuousGeomCurves = new QCheckBox(verticalLayoutWidget);
        checkBox_DisplayContinuousGeomCurves->setObjectName(QString::fromUtf8("checkBox_DisplayContinuousGeomCurves"));

        verticalLayout->addWidget(checkBox_DisplayContinuousGeomCurves);

        checkBox_DisplayDiscretizedGeomCurves = new QCheckBox(verticalLayoutWidget);
        checkBox_DisplayDiscretizedGeomCurves->setObjectName(QString::fromUtf8("checkBox_DisplayDiscretizedGeomCurves"));

        verticalLayout->addWidget(checkBox_DisplayDiscretizedGeomCurves);

        label_3 = new QLabel(frame);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(70, 130, 121, 16));
        spinBox_CurrentEdgeID = new QSpinBox(frame);
        spinBox_CurrentEdgeID->setObjectName(QString::fromUtf8("spinBox_CurrentEdgeID"));
        spinBox_CurrentEdgeID->setGeometry(QRect(140, 310, 111, 31));
        label_8 = new QLabel(frame);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(10, 360, 111, 16));
        lineEdit_NumOfEdgeSegments = new QLineEdit(frame);
        lineEdit_NumOfEdgeSegments->setObjectName(QString::fromUtf8("lineEdit_NumOfEdgeSegments"));
        lineEdit_NumOfEdgeSegments->setGeometry(QRect(140, 350, 111, 31));
        label_7 = new QLabel(frame);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(10, 320, 121, 16));
        gridLayoutWidget_2 = new QWidget(frame);
        gridLayoutWidget_2->setObjectName(QString::fromUtf8("gridLayoutWidget_2"));
        gridLayoutWidget_2->setGeometry(QRect(10, 150, 241, 151));
        gridLayout_2 = new QGridLayout(gridLayoutWidget_2);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setContentsMargins(0, 0, 0, 0);
        lineEdit_MinEdgeLength = new QLineEdit(gridLayoutWidget_2);
        lineEdit_MinEdgeLength->setObjectName(QString::fromUtf8("lineEdit_MinEdgeLength"));

        gridLayout_2->addWidget(lineEdit_MinEdgeLength, 0, 2, 1, 1);

        label_5 = new QLabel(gridLayoutWidget_2);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_2->addWidget(label_5, 1, 0, 1, 1);

        lineEdit_MaxEdgeLength = new QLineEdit(gridLayoutWidget_2);
        lineEdit_MaxEdgeLength->setObjectName(QString::fromUtf8("lineEdit_MaxEdgeLength"));

        gridLayout_2->addWidget(lineEdit_MaxEdgeLength, 1, 2, 1, 1);

        label_4 = new QLabel(gridLayoutWidget_2);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_2->addWidget(label_4, 0, 0, 1, 1);

        label_6 = new QLabel(gridLayoutWidget_2);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_2->addWidget(label_6, 2, 0, 1, 1);

        lineEdit_MeanEdgeLength = new QLineEdit(gridLayoutWidget_2);
        lineEdit_MeanEdgeLength->setObjectName(QString::fromUtf8("lineEdit_MeanEdgeLength"));

        gridLayout_2->addWidget(lineEdit_MeanEdgeLength, 2, 2, 1, 1);

        line_2 = new QFrame(frame);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setGeometry(QRect(0, 410, 311, 16));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        label = new QLabel(frame);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(80, 430, 111, 21));
        radioButton_DisplayUVSurface = new QRadioButton(frame);
        radioButton_DisplayUVSurface->setObjectName(QString::fromUtf8("radioButton_DisplayUVSurface"));
        radioButton_DisplayUVSurface->setGeometry(QRect(10, 500, 106, 23));
        radioButton_DisplayXYZSurface = new QRadioButton(frame);
        radioButton_DisplayXYZSurface->setObjectName(QString::fromUtf8("radioButton_DisplayXYZSurface"));
        radioButton_DisplayXYZSurface->setGeometry(QRect(140, 500, 106, 23));
        label_2 = new QLabel(frame);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 470, 111, 16));
        spinBox_CurrentFaceID = new QSpinBox(frame);
        spinBox_CurrentFaceID->setObjectName(QString::fromUtf8("spinBox_CurrentFaceID"));
        spinBox_CurrentFaceID->setGeometry(QRect(140, 460, 111, 31));
        pushButton_GenerateSurface = new QPushButton(frame);
        pushButton_GenerateSurface->setObjectName(QString::fromUtf8("pushButton_GenerateSurface"));
        pushButton_GenerateSurface->setGeometry(QRect(10, 580, 81, 31));
        pushButton_AcceptSurface = new QPushButton(frame);
        pushButton_AcceptSurface->setObjectName(QString::fromUtf8("pushButton_AcceptSurface"));
        pushButton_AcceptSurface->setGeometry(QRect(100, 580, 71, 31));
        pushButton_RejectSurface = new QPushButton(frame);
        pushButton_RejectSurface->setObjectName(QString::fromUtf8("pushButton_RejectSurface"));
        pushButton_RejectSurface->setGeometry(QRect(190, 580, 61, 31));
        label_9 = new QLabel(frame);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(70, 10, 121, 21));
        checkBox_displayCurrentEdge = new QCheckBox(frame);
        checkBox_displayCurrentEdge->setObjectName(QString::fromUtf8("checkBox_displayCurrentEdge"));
        checkBox_displayCurrentEdge->setGeometry(QRect(10, 390, 231, 21));
        checkBox_DisplayCurrentFace = new QCheckBox(frame);
        checkBox_DisplayCurrentFace->setObjectName(QString::fromUtf8("checkBox_DisplayCurrentFace"));
        checkBox_DisplayCurrentFace->setGeometry(QRect(10, 540, 221, 21));
        line = new QFrame(frame);
        line->setObjectName(QString::fromUtf8("line"));
        line->setGeometry(QRect(20, 130, 229, 3));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1151, 25));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        menubar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(actionExit);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        checkBox_displayGeometricNodes->setText(QApplication::translate("MainWindow", "Display Nodes", 0, QApplication::UnicodeUTF8));
        checkBox_DisplayContinuousGeomCurves->setText(QApplication::translate("MainWindow", "Display Continuous Curves", 0, QApplication::UnicodeUTF8));
        checkBox_DisplayDiscretizedGeomCurves->setText(QApplication::translate("MainWindow", "Display Discretized Curves", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Discretize Edge", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "Num Segments", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("MainWindow", "Current Edge ID", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Max Edge Length ", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Min Edge Length", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "Mean Edge Length", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Discretize Face", 0, QApplication::UnicodeUTF8));
        radioButton_DisplayUVSurface->setText(QApplication::translate("MainWindow", "Display UV", 0, QApplication::UnicodeUTF8));
        radioButton_DisplayXYZSurface->setText(QApplication::translate("MainWindow", "Display XYZ", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Current Face ID", 0, QApplication::UnicodeUTF8));
        pushButton_GenerateSurface->setText(QApplication::translate("MainWindow", "Generate", 0, QApplication::UnicodeUTF8));
        pushButton_AcceptSurface->setText(QApplication::translate("MainWindow", "Accept", 0, QApplication::UnicodeUTF8));
        pushButton_RejectSurface->setText(QApplication::translate("MainWindow", "Reject", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("MainWindow", "Geometric Model", 0, QApplication::UnicodeUTF8));
        checkBox_displayCurrentEdge->setText(QApplication::translate("MainWindow", "Display Current Edge", 0, QApplication::UnicodeUTF8));
        checkBox_DisplayCurrentFace->setText(QApplication::translate("MainWindow", "Display Current Face", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ITAPS_H
