#include <iostream>
#include <Windows.h>
#include <stdio.h>
#include <stdlib.h>

#include <opencv\cv.h>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/nonfree/nonfree.hpp>

#include <cmath>

using namespace cv;
using namespace std;
Mat C, Cinv, CROI, CirROI, CTriROI, Tinv, H;
Mat Window=(Mat_<uchar>(3,3)<<0,1,0,1,1,1,0,1,0);

int ho,ve; double prev_area=0.0;
int maxLength, maxLengthROI;
int count_area=0, count_ave=0;

int max_thresh=240;//===> was 200, changed for experiment on 5 oct

double average_area=0.0, sum_area=0.0, average_width=0.0, sum_width=0.0, average_height=0.0, sum_height=0.0;
Point2f src[4], dst[4], srcROI[4], dstROI[4];
Point2f MoveTLpt,CMoveTLpt;
Point2f proPt, cor;
Point2f rectCentre, prev_rectCentre, rectAngle, prev_rectAngle, rectPoint, prev_rectPoint;

double 	rectSize_sumW=0.0,rectSize_sumH=0.0, average_W=0.0,average_H=0.0;

bool call, method;
bool Method;
bool circular, triangular;

void drawLine(Vec2f line, Mat &img, Scalar rgb=CV_RGB(0,0,255));
void mergeRelatedLinesCS(vector<Vec2f>*lines, Mat &img);
void mergeRelatedLinesPC(vector<Vec2f>*lines, Mat &img);
void mergeRelatedLinesROI(vector<Vec2f>*lines, Mat &img);

void drawlines(vector<Vec2f> &lines,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge);
void drawlinesROI(vector<Vec2f> &lines,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge);

void intersectionPoints(Mat &Border,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge, Point2f &ptTopLeft, Point2f &ptTopRight, Point2f &ptBotLeft, Point2f &ptBotRight);
int findMaxLength(int &maxLength, Point2f &ptTopLeft, Point2f &ptTopRight, Point2f &ptBotLeft, Point2f &ptBotRight);
void floodFilling(Mat &Img);

void CameraScreen (Mat& Screen, Mat&C);
void ProjectorCamera (Mat &Dst, Mat& Src, Point2f srcpt[4], Point2f dstpt[4]);
void GetROIinScreen(Mat &Screen, Mat & CROI, Mat C);

double FindMax(double a, double b, double c, double d);
double FindMin(double a,double b,double c, double d);
bool SameSign(int x, int y);
void GetDesktopResolution(int& horizontal, int& vertical);
void SetExposure(VideoCapture Cap);

int maxH=110, minH=0, H_slide; double Hu=110;
int maxS=121, minS=0, S_slide; double Sa=121;
int maxV=231, minV=0, V_slide; double Va=231;

void on_H(int, void*);
void on_S(int, void*);
void on_V(int, void*);
double rectCentreDiff(Point2f cen, Point2f pcen);
double rectArea(Point2f s[4]);
int distanceCheck(Point2f s[4]);
void ShowSrcImage(Mat Src);
void equaliseGray(Mat & Gray);
void sortCorners(Point2f r[4], int img_w, int img_h);
void allocateCorners(Point2f rectpoints[4], Mat Rundistorted);
void identifyShape(Mat Screen, bool& circular, bool& triangular);


Mat getTransform2(Point2f srcpt[16], Point2f dstpt[16]);
Mat getTransform(Point2f srcpt[4], Point2f dstpt[4]);
Mat HSV(VideoCapture& Cap, Mat& Input);
Mat DetectBorder(Mat Rundistorted);
Mat DetectCircle(Mat Crop);
Mat DetectTriangle(Mat Crop);
Mat translateImg(Mat &img, int offsetx, int offsety);
Mat MakeProjection (Mat &Screen, Mat T,Mat C, bool Method);
Mat WarpImage(Mat Image, Mat P, double ScaleX, double ScaleY, bool circular,bool triangular, Mat T);
Mat MakeCircleProjection(Mat Screen, Mat C);
Mat MakeTriProjection(Mat Screen, Mat C);
Mat Make_n_Warp(Mat Screen, Mat T, Mat C, Mat pic);
Mat Tri_Make_n_Warp(Mat Screen, Mat T, Mat C, Mat pic);
Mat Capturing(Mat Src);
Mat ProCamMatrix(bool debug);

void filterNoise(Point2f s[4]);
void filterInit(Point2f s[4]);

void RfilterNoise(Point2f s[4]);
void RfilterInit(Point2f s[4]);

void CenfilterNoise(Point2f &s);
void CenfilterInit(Point2f &s);

void AngfilterNoise(Point2f &s);
void AngfilterInit(Point2f &s);

KalmanFilter KF1(4,2,0);
KalmanFilter KF2(4,2,0);
KalmanFilter KF3(4,2,0);
KalmanFilter KF4(4,2,0);

Mat_<float> measure1(2,1);
Mat_<float> measure2(2,1);
Mat_<float> measure3(2,1);
Mat_<float> measure4(2,1);

KalmanFilter CenKF(6,2,0);
Mat_<float> Cenmeasure(2,1);

KalmanFilter AngKF(6,2,0);
Mat_<float> Angmeasure(2,1);

KalmanFilter RKF1(6,2,0);
KalmanFilter RKF2(6,2,0);
KalmanFilter RKF3(6,2,0);
KalmanFilter RKF4(6,2,0);

Mat_<float> Rmeasure1(2,1);
Mat_<float> Rmeasure2(2,1);
Mat_<float> Rmeasure3(2,1);
Mat_<float> Rmeasure4(2,1);
