#include "header.h"

void equaliseGray(Mat & Gray)
{
	double alpha=2.2;	//1.0-3.0
	int beta=50;		//0-100
	Scalar mymean, stddev;
	double stdv;

	meanStdDev(Gray, mymean, stddev);
	stdv= stddev.val[0];
	if(stdv<20)
	{
		Gray.convertTo(Gray, -1, alpha, beta);
	}
}

void allocateCorners(Point2f rectpoints[4], Mat Rundistorted)
{
	Point2f ptTopLeftROI, ptTopRightROI, ptBotLeftROI, ptBotRightROI;
	sortCorners(rectpoints, Rundistorted.cols, Rundistorted.rows);

	ptTopLeftROI=rectpoints[0];	ptTopRightROI=rectpoints[1];
	ptBotLeftROI=rectpoints[3];	ptBotRightROI=rectpoints[2];
  //cout<<"ptTopLeftROI= ["<<ptTopLeftROI.x<<","<<ptTopLeftROI.y<<"]"<<endl;
	//cout<<"ptTopRightROI= ["<<ptTopRightROI.x<<","<<ptTopRightROI.y<<"]"<<endl;
	//cout<<"ptBotRightROI= ["<<ptBotRightROI.x<<","<<ptBotRightROI.y<<"]"<<endl;
	//cout<<"ptBotLeftROI= ["<<ptBotLeftROI.x<<","<<ptBotLeftROI.y<<"]"<<endl;


	/*ptBotRightROI=rectpoints[0];	ptBotLeftROI=rectpoints[1];
	ptTopRightROI=rectpoints[3];	ptTopLeftROI=rectpoints[2];*///===> wr bt updw
	srcROI[0]=ptTopLeftROI;		srcROI[1]=ptTopRightROI; 
	srcROI[2]=ptBotRightROI;	srcROI[3]=ptBotLeftROI;
	findMaxLength(maxLengthROI,ptTopLeftROI,ptTopRightROI,ptBotLeftROI,ptBotRightROI);
	//(maxlength: maxlength*0.70707), where aspect ratio is 4:3, for width: height
	
	dstROI[0]=Point2f(0,0);		dstROI[1]=Point2f(maxLengthROI-1,0);
	dstROI[2]=Point2f(maxLengthROI-1,floor(maxLengthROI*0.70707)-1);
	dstROI[3]=Point2f(0,floor(maxLengthROI*0.70707)-1);
}

void sortCorners(Point2f r[4], int img_w, int img_h)
{
	Point2f sp_l[4], sp_r[4];
	sp_l[0]=r[0];	sp_l[1]=r[1];
	sp_l[2]=r[2];	sp_l[3]=r[3];
	
	double min_dist=10000, max_dist=0;
	double min_dist_r=10000, max_dist_r=0;
	for (int i =0;i<4;i++)
	{
		double dist= sqrt((sp_l[i].x*sp_l[i].x) + (sp_l[i].y*sp_l[i].y));
		//
		if(dist<min_dist)
		{
			min_dist=dist;
			r[0]=sp_l[i];
		}
		/*if(dist>max_dist)
		{
			max_dist=dist;
			r[2]=sp_l[i];
		}
		if(dist_r<min_dist_r)
		{
			min_dist_r=dist_r;
			r[3]=sp_l[i];
		}
		if(dist_r>max_dist_r)
		{
			max_dist_r=dist_r;
			r[1]=sp_l[i];
		}*/
	}
	double min_ang=90, max_ang=0;
	for (int i =0;i<4;i++)
	{
		if(sp_l[i].y==r[0].y && sp_l[i].x==r[0].x)
			continue;
		float angle= atan2(sp_l[i].y,sp_l[i].x)*180/3.14159265;
		if(angle<=min_ang)
		{
			min_ang=angle;
			r[1]=sp_l[i];
		}
		if(angle>=max_ang)
		{
			max_ang=angle;
			r[3]=sp_l[i];
		}
	}
	for (int i =0;i<4;i++)
	{
		if(sp_l[i]!=r[0] && sp_l[i]!=r[1] && sp_l[i]!=r[3])
		{
			r[2]=sp_l[i];
		}
	}
}

Mat translateImg(Mat &img, int offsetx, int offsety)
{
    Mat trans_mat = (Mat_<double>(2,3) << 1, 0, offsetx, 0, 1, offsety);
    warpAffine(img,img,trans_mat,img.size());
    return trans_mat;
}

int distanceCheck(Point2f s[4])
{
	double a,b,c,d;
	int fa;

	a= sqrt((s[0].x-s[1].x)*(s[0].x-s[1].x) + (s[0].y-s[1].y)*(s[0].y-s[1].y));
	b= sqrt((s[1].x-s[2].x)*(s[1].x-s[2].x) + (s[1].y-s[2].y)*(s[1].y-s[2].y));
	c= sqrt((s[2].x-s[3].x)*(s[2].x-s[3].x) + (s[2].y-s[3].y)*(s[2].y-s[3].y));
	d= sqrt((s[3].x-s[0].x)*(s[3].x-s[0].x) + (s[3].y-s[0].y)*(s[3].y-s[0].y));
	if (a>10 && b>10 && c>10 && d>10)
		fa=1;
	else
		fa=0;

	return fa;
}

double rectArea(Point2f s[4])
{
	Point2f r[4];
	r[0]=s[0]; r[1]=s[1]; r[2]=s[2]; r[3]=s[3];
	double rect_area= sqrt((r[1].x-r[0].x)*(r[1].x-r[0].x)+(r[1].y-r[0].y)*(r[1].y-r[0].y))*
		sqrt((r[3].x-r[0].x)*(r[3].x-r[0].x)+(r[3].y-r[0].y)*(r[3].y-r[0].y));
	return rect_area;
}

double rectCentreDiff(Point2f cen, Point2f pcen)
{
	double diff=sqrt((cen.x-pcen.x)*(cen.x-pcen.x)+(cen.y-pcen.y)*(cen.y-pcen.y));
	return diff;
}

double rectAngleDiff(Point2f ang, Point2f pang)
{
	double diff=sqrt((ang.x-pang.x)*(ang.x-pang.x)+(ang.y-pang.y)*(ang.y-pang.y));
	return diff;
}

void drawLine(Vec2f line, Mat &img, Scalar rgb)
{
 if(line[1]!=0)
 {
      float m=-1/tan(line[1]);
      float c=line[0]/sin(line[1]);
      cv::line(img,Point(0,c),Point(img.size().width,m*img.size().width+c),rgb);
 }
 else
 {
    cv::line(img,Point(line[0],0),Point(line[0],img.size().height),rgb);
 }
}

void mergeRelatedLinesCS(vector<Vec2f>*lines, Mat &img)
{
 vector<Vec2f>::iterator current;
 for( current=lines->begin(); current!=lines->end();current++)
 {
  if((*current)[0]==0 && (*current)[1]==-100)
   continue;
  float p1=(*current)[0];
  float theta1=(*current)[1];
  Point pt1current, pt2current;
  if(theta1>CV_PI*45/180 && theta1<CV_PI*135/180)
  {
       pt1current.x=0;
       pt1current.y=p1/sin(theta1);
       pt2current.x=img.size().width;
       pt2current.y=-pt2current.x/tan(theta1) + p1/sin(theta1);
  }
  else
  {
       pt1current.y=0;
       pt1current.x=p1/cos(theta1);
       pt2current.y=img.size().height;
       pt2current.x=-pt2current.y/tan(theta1) + p1/cos(theta1);
  }

  vector<Vec2f>::iterator pos;
  for (pos=lines->begin();pos!=lines->end();pos++)
  {
      if(*current==*pos)
        continue;
      if(fabs((*pos)[0]-(*current)[0])<20 && fabs((*pos)[1]-(*current)[1])<CV_PI*10/180)
      {
            float p=(*pos)[0];
            float theta=(*pos)[1];
            Point pt1, pt2;
            if((*pos)[1]>CV_PI*45/180 && (*pos)[1]<CV_PI*135/180)
            {
                 pt1.x=0;
                 pt1.y=p/sin(theta);
                 pt2.x=img.size().width;
                 pt2.y=-pt2.x/tan(theta) + p/sin(theta);
            }
            else 
            {
                 pt1.y=0;
                 pt1.x=p/cos(theta);
                 pt2.y=img.size().height;
                 pt2.x=-pt2.y/tan(theta) + p/cos(theta);
            }
            if(((double)(pt1.x-pt1current.x)*(pt1.x-pt1current.x) + (pt1.y-pt1current.y)*(pt1.y-pt1current.y)<128*128) &&
             ((double)(pt2.x-pt2current.x)*(pt2.x-pt2current.x) +(pt2.y-pt2current.y)*(pt2.y-pt2current.y)<128*128))
            {
                 (*current)[0]=((*current)[0]+(*pos)[0])/2;
                 (*current)[1]=((*current)[1]+(*pos)[1])/2;
                 (*pos)[0]=0;
                 (*pos)[1]=-100;
            }
        }
      }
   }
}

void mergeRelatedLinesROI(vector<Vec2f>*lines, Mat &img)
{
 vector<Vec2f>::iterator current;
 for( current=lines->begin(); current!=lines->end();current++)
 {
  if((*current)[0]==0 && (*current)[1]==-100)
   continue;
  float p1=(*current)[0];
  float theta1=(*current)[1];
  Point pt1current, pt2current;
  if(theta1>CV_PI*45/180 && theta1<CV_PI*135/180)
  {
       pt1current.x=0;
       pt1current.y=p1/sin(theta1);
       pt2current.x=img.size().width;
       pt2current.y=-pt2current.x/tan(theta1) + p1/sin(theta1);
  }
  else
  {
       pt1current.y=0;
       pt1current.x=p1/cos(theta1);
       pt2current.y=img.size().height;
       pt2current.x=-pt2current.y/tan(theta1) + p1/cos(theta1);
  }

  vector<Vec2f>::iterator pos;
  for (pos=lines->begin();pos!=lines->end();pos++)
  {
      if(*current==*pos)
        continue;
      if(fabs((*pos)[0]-(*current)[0])<30 && fabs((*pos)[1]-(*current)[1])<CV_PI*30/180)
      {
            float p=(*pos)[0];
            float theta=(*pos)[1];
            Point pt1, pt2;
            if((*pos)[1]>CV_PI*45/180 && (*pos)[1]<CV_PI*135/180)
            {
                 pt1.x=0;
                 pt1.y=p/sin(theta);
                 pt2.x=img.size().width;
                 pt2.y=-pt2.x/tan(theta) + p/sin(theta);
            }
            else 
            {
                 pt1.y=0;
                 pt1.x=p/cos(theta);
                 pt2.y=img.size().height;
                 pt2.x=-pt2.y/tan(theta) + p/cos(theta);
            }
            if(((double)(pt1.x-pt1current.x)*(pt1.x-pt1current.x) + (pt1.y-pt1current.y)*(pt1.y-pt1current.y)<128*128) &&
             ((double)(pt2.x-pt2current.x)*(pt2.x-pt2current.x) +(pt2.y-pt2current.y)*(pt2.y-pt2current.y)<128*128))
            {
                 (*current)[0]=((*current)[0]+(*pos)[0])/2;
                 (*current)[1]=((*current)[1]+(*pos)[1])/2;
                 (*pos)[0]=0;
                 (*pos)[1]=-100;
            }
        }
      }
    }
}

void mergeRelatedLinesPC(vector<Vec2f>*lines, Mat &img)
{
 vector<Vec2f>::iterator current;
 for( current=lines->begin(); current!=lines->end();current++)
 {
  if((*current)[0]==0 && (*current)[1]==-100)
   continue;
  float p1=(*current)[0];
  float theta1=(*current)[1];
  Point pt1current, pt2current;
  if(theta1>CV_PI*45/180 && theta1<CV_PI*135/180)
  {
       pt1current.x=0;
       pt1current.y=p1/sin(theta1);
       pt2current.x=img.size().width;
       pt2current.y=-pt2current.x/tan(theta1) + p1/sin(theta1);
  }
  else
  {
       pt1current.y=0;
       pt1current.x=p1/cos(theta1);
       pt2current.y=img.size().height;
       pt2current.x=-pt2current.y/tan(theta1) + p1/cos(theta1);
  }

  vector<Vec2f>::iterator pos;
  for (pos=lines->begin();pos!=lines->end();pos++)
  {
      if(*current==*pos)
        continue;
      if(fabs((*pos)[0]-(*current)[0])<100 && fabs((*pos)[1]-(*current)[1])<CV_PI*30/180)
      {
            float p=(*pos)[0];
            float theta=(*pos)[1];
            Point pt1, pt2;
            if((*pos)[1]>CV_PI*45/180 && (*pos)[1]<CV_PI*135/180)
            {
                 pt1.x=0;
                 pt1.y=p/sin(theta);
                 pt2.x=img.size().width;
                 pt2.y=-pt2.x/tan(theta) + p/sin(theta);
            }
            else 
            {
                 pt1.y=0;
                 pt1.x=p/cos(theta);
                 pt2.y=img.size().height;
                 pt2.x=-pt2.y/tan(theta) + p/cos(theta);
            }
            if(((double)(pt1.x-pt1current.x)*(pt1.x-pt1current.x) + (pt1.y-pt1current.y)*(pt1.y-pt1current.y)<150*150) &&
             ((double)(pt2.x-pt2current.x)*(pt2.x-pt2current.x) +(pt2.y-pt2current.y)*(pt2.y-pt2current.y)<150*150))
            {
                 (*current)[0]=((*current)[0]+(*pos)[0])/2;
                 (*current)[1]=((*current)[1]+(*pos)[1])/2;
                 (*pos)[0]=0;
                 (*pos)[1]=-100;
            }
        }
      }
    }
}

void drawlines(vector<Vec2f> &lines,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge)
{
    double topYintercept=100000, topXintercept=0; double bottomYintercept=0, bottomXintercept=0;
    double leftXintercept=100000, leftYintercept=0; double rightXintercept=0, rightYintercept=0;
    for(int i=0; i<lines.size(); i++)
    {
        Vec2f current=lines[i];
        float p=current[0];
        float theta=current[1];
        if(p==0 && theta==-100)
            continue;
        double xintercept, yintercept;
        xintercept= p/cos(theta);
        yintercept= p/(cos(theta)*sin(theta));
        if(theta>CV_PI*70/180 && theta<CV_PI*100/180)
        {
            if(p<topEdge[0])
                topEdge=current;
            if(p>bottomEdge[0])
                bottomEdge=current;
        }
        else if (theta<CV_PI*10/180 || theta>CV_PI*170/180)
        {
            if(xintercept>rightXintercept)
            {
                rightEdge=current;
                rightXintercept=xintercept;
            }
            else if(xintercept<=leftXintercept)
            {
                leftEdge=current;
                leftXintercept=xintercept;
            }
        }
    }
}

void drawlinesROI(vector<Vec2f> &lines,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge)
{
    double topYintercept=100000, topXintercept=0; double bottomYintercept=0, bottomXintercept=0;
    double leftXintercept=100000, leftYintercept=0; double rightXintercept=0, rightYintercept=0;
	int count =0;
    for(int i=0; i<lines.size(); i++)
    {
        Vec2f current=lines[i];
        float p=current[0];
        float theta=current[1];
        if(p==0 && theta==-100)
            continue;
        double xintercept, yintercept;
        xintercept= p/cos(theta);
        yintercept= p/(cos(theta)*sin(theta));
        if(theta>CV_PI*45/180 && theta<CV_PI*135/180)
        {
            if(p<topEdge[0])
                topEdge=current;
            if(p>bottomEdge[0])
                bottomEdge=current;
        }
        else if (theta<CV_PI*45/180 || theta>CV_PI*135/180)
        {
            if(xintercept>rightXintercept)
            {
				if(rightXintercept!=0 && abs(rightXintercept-xintercept)>=200)
				{
					leftEdge=rightEdge;
					leftXintercept=rightXintercept;
					rightEdge=current;
					rightXintercept=xintercept;
				}
				else
				{
					rightEdge=current;
					rightXintercept=xintercept;
				}
            }
            else if(xintercept<=leftXintercept)
            {
                leftEdge=current;
                leftXintercept=xintercept;
            }
        }
    }
}

void intersectionPoints(Mat &Border,Vec2f &topEdge, Vec2f &bottomEdge, Vec2f &rightEdge, Vec2f &leftEdge, Point2f &ptTopLeft, Point2f &ptTopRight, Point2f &ptBotLeft, Point2f &ptBotRight)
{
     Point left1,left2, right1,right2, bot1,bot2, top1,top2;
     int height=Border.size().height;
     int width=Border.size().width;

     if(leftEdge[1]!=0)
     {
          left1.x=0;  left1.y=leftEdge[0]/sin(leftEdge[1]);
          left2.x=width; left2.y=-left2.x/tan(leftEdge[1])+left1.y;
     }
     else
     {
          left1.y=0;   left1.x=leftEdge[0]/cos(leftEdge[1]);
          left2.y=height; left2.x=left1.x-height*tan(leftEdge[1]);
     }
     if(rightEdge[1]!=0)
     {
          right1.x=0;     right1.y=rightEdge[0]/sin(rightEdge[1]);
          right2.x=width; right2.y=-right2.x/tan(rightEdge[1])+right1.y;
     }
     else
     {
          right1.y=0;      right1.x=rightEdge[0]/cos(rightEdge[1]);
          right2.y=height; right2.x=right1.x-height*tan(rightEdge[1]);
     }
     bot1.x=0; bot1.y=bottomEdge[0]/sin(bottomEdge[1]);
     bot2.x=width; bot2.y=-bot2.x/tan(bottomEdge[1])+bot1.y;

     top1.x=0; top1.y=topEdge[0]/sin(topEdge[1]);
     top2.x=width; top2.y=-top2.x/tan(topEdge[1])+top1.y;

     double leftA=left2.y-left1.y; double leftB=left1.x-left2.x;
     double leftC=leftA*left1.x+leftB*left1.y;
 
     double rightA=right2.y-right1.y; double rightB=right1.x-right2.x;
     double rightC=rightA*right1.x+rightB*right1.y;

     double topA=top2.y-top1.y; double topB=top1.x-top2.x;
     double topC=topA*top1.x+topB*top1.y;

     double botA=bot2.y-bot1.y; double botB=bot1.x-bot2.x;
     double botC=botA*bot1.x+botB*bot1.y;

     double detTopLeft=leftA*topB-leftB*topA;
     ptTopLeft=cvPoint((topB*leftC-leftB*topC)/detTopLeft,(leftA*topC-topA*leftC)/detTopLeft);
     double detTopRight=rightA*topB-rightB*topA;
     ptTopRight=cvPoint((topB*rightC-rightB*topC)/detTopRight,(rightA*topC-topA*rightC)/detTopRight);
 
     double detBotRight=rightA*botB-rightB*botA;
     ptBotRight=cvPoint((botB*rightC-rightB*botC)/detBotRight,(rightA*botC-botA*rightC)/detBotRight);
     double detBotLeft=leftA*botB-leftB*botA;
     ptBotLeft=cvPoint((botB*leftC-leftB*botC)/detBotLeft,(leftA*botC-botA*leftC)/detBotLeft);
}

int findMaxLength(int &maxLength, Point2f &ptTopLeft, Point2f &ptTopRight, Point2f &ptBotLeft, Point2f &ptBotRight)
{
     maxLength =(ptBotLeft.x-ptBotRight.x)*(ptBotLeft.x-ptBotRight.x)+(ptBotLeft.y-ptBotRight.y)*(ptBotLeft.y-ptBotRight.y);
     int temp= (ptTopRight.x-ptBotRight.x)*(ptTopRight.x-ptBotRight.x)+(ptTopRight.y-ptBotRight.y)*(ptTopRight.y-ptBotRight.y);
     if(temp>maxLength)
      maxLength=temp;
     temp=(ptTopRight.x-ptTopLeft.x)*(ptTopRight.x-ptTopLeft.x)+(ptTopRight.y-ptTopLeft.y)*(ptTopRight.y-ptTopLeft.y);
     if(temp>maxLength)
      maxLength=temp;
     temp=(ptBotLeft.x-ptTopLeft.x)*(ptBotLeft.x-ptTopLeft.x)+(ptBotLeft.y-ptTopLeft.y)*(ptBotLeft.y-ptTopLeft.y);
     if(temp>maxLength)
      maxLength=temp;
     return maxLength=sqrt((double)maxLength);
}

void floodFilling(Mat &Img)
{
	int max=-1;
	Point MaxPt;
	for ( int y=0; y< Img.size().height;y++)
    {
        uchar *row=Img.ptr(y);
        for( int x=0;x<Img.size().width;x++)
        {
           if (row[x]>128)
           {
             int area=floodFill(Img,Point(x,y),CV_RGB(0,0,64));
             if(area>max)
             {
               MaxPt=Point(x,y);
               max=area;
             }
           }
        }
    }
    floodFill(Img, MaxPt, CV_RGB(255,255,255));
    for(int y=0;y<Img.size().height;y++)
    {
         uchar * row=Img.ptr(y);
         for(int x=0;x<Img.size().width;x++)
         {
             if(row[x]==64 && x!=MaxPt.x && y!=MaxPt.y)
             {
                 int area=floodFill(Img,Point(x,y),CV_RGB(0,0,0));
             }
         }
    }
}

double FindMax(double a, double b, double c, double d)
 {
	 double Max1, Max2, Max3, Max4, Max;
	 Max1 = max(a,b);
	 Max2 = max(c,d);
	 Max3 = max(b,c);
	 Max4 = max(Max1,Max3);
	 Max =  max(Max2,Max4);
	 return Max;
 }

double FindMin(double a,double b,double c, double d)
 {
	 double Min1, Min2, Min3, Min4, Min;
	 Min1 = min(a,b);
	 Min2 = min(c,d);
	 Min3 = min(b,c);
	 Min4 = min(Min1,Min3);
	 Min  = min(Min2,Min4);
	 return Min;
 }

bool SameSign(int x, int y)
{
    return (x >= 0) ^ (y < 0);
}

void GetDesktopResolution(int& horizontal, int& vertical)
{
   RECT desktop;
   const HWND hDesktop = GetDesktopWindow();
   GetWindowRect(hDesktop, &desktop);
   horizontal = desktop.right;
   vertical = desktop.bottom;
}

void CameraScreen (Mat& Screen, Mat &C)
{
	Mat Border, Gray, Undistorted, Window;
	vector<Vec2f> linesCS;
	cvtColor(Screen,Gray,CV_BGR2GRAY);
	Vec2f topEdgeCS, bottomEdgeCS, leftEdgeCS, rightEdgeCS; 
	Point2f ptTopLeftCS, ptTopRightCS, ptBotLeftCS, ptBotRightCS;
	double alpha=2.2;	//1.0-3.0
	int beta=50;		//0-100
	Scalar mymean, stddev;
	double stdv;

	meanStdDev(Gray, mymean, stddev);
	stdv= stddev.val[0];
	if(stdv<20)
	{
		Gray.convertTo(Gray, -1, alpha, beta);
	}
	//imshow("gray",Gray);
	//waitKey(1);
    GaussianBlur(Gray, Gray, Size(5,5),0);
    adaptiveThreshold(Gray,Border,max_thresh,ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,11,6);
	if(GetAsyncKeyState(VK_DOWN))
		max_thresh=max_thresh--;
	if(GetAsyncKeyState(VK_UP))
		max_thresh=max_thresh++;

	//threshold(Gray,Border,20,255,THRESH_TOZERO_INV);
	//imshow("Border1",Border);
	//waitKey(1);
	bitwise_not(Border,Border);


    dilate(Border,Border,Window);
	floodFilling(Border);
    erode(Border,Border,Window);
	//imshow("Border",Border);
	//waitKey(1);

	//Mat Border2= Mat::zeros(Border.rows, Border.cols, Border.type());
	//Scalar color1 = Scalar(255,255,255);
	//Border.copyTo(Border2);

	Canny(Border,Border,20,250,3);
	imshow("Border",Border);
	waitKey(1);

	//vector<vector<Point>> contourst;
	//vector<Vec4i> hiert;
	//int largestareat=0; int largestindext=0;
	//findContours(Border2,contourst,hiert,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE
	//	
	//vector<vector<Point>>hullt( contourst.size() );
	//					 
	//for( int i = 0; i < contourst.size(); i++ )
	//convexHull(Mat(contourst[i]), hullt[i],false);
	//
	//for( int i = 0; i < hullt.size(); i++ )
	//{ 
	//	double iter = contourArea(hullt[i],false);
	//	if (iter>largestareat)
	//	{
	//		largestareat=iter;
	//		largestindext=i;
	//	}	 					
	//}

	//drawContours( Border2, hullt, 0, color1, 1, 8, vector<Vec4i>(), 0, Point() );
	//imshow("Border2",Border2);
	//waitKey(1);

	//Border2.copyTo(Border);
       
	HoughLines(Border,linesCS,1,CV_PI/180,50);//=============> 75 changes to 50
    mergeRelatedLinesCS(&linesCS,Gray);
    topEdgeCS=Vec2f(1000,1000); bottomEdgeCS=Vec2f(-1000,-1000);
    leftEdgeCS=Vec2f(1000,1000); rightEdgeCS=Vec2f(-1000,-1000);

	drawlines(linesCS,topEdgeCS,bottomEdgeCS,rightEdgeCS,leftEdgeCS);
    intersectionPoints(Border,topEdgeCS,bottomEdgeCS,rightEdgeCS,leftEdgeCS,ptTopLeftCS,ptTopRightCS,ptBotLeftCS,ptBotRightCS);
	//-------------------------------------------//
	//drawLine(topEdgeCS, Screen, Scalar(255,255,255));
	//drawLine(bottomEdgeCS, Screen, Scalar(255,255,255));
	//drawLine(rightEdgeCS, Screen, Scalar(255,255,255));
	//drawLine(leftEdgeCS, Screen, Scalar(255,255,255));
	//imshow("Screen",Screen);
	//waitKey(1);
	//-------------------------------------------//

	findMaxLength(maxLength,ptTopLeftCS,ptTopRightCS,ptBotLeftCS,ptBotRightCS);
    //(maxlength: maxlength*3/4), where aspect ratio is 4:3, for width: height
    src[0]=ptTopLeftCS; dst[0]=Point2f(0,0);
    src[1]=ptTopRightCS; dst[1]=Point2f(maxLength-1,0);
    src[2]=ptBotRightCS; dst[2]=Point2f(maxLength-1,floor(maxLength*0.75)-1);
    src[3]=ptBotLeftCS; dst[3]=Point2f(0,floor(maxLength*0.75)-1);

	filterNoise(src);
	CMoveTLpt=src[0];
 
	C=getPerspectiveTransform(src,dst);
	Cinv = getPerspectiveTransform(dst,src);

	//Undistorted = Mat(Size(maxLength,floor(maxLength*0.75)),CV_8UC1);
    cv::warpPerspective(Gray,Undistorted,cv::getPerspectiveTransform(src,dst),Size(maxLength,floor(maxLength*0.75)));
    //cout<<"The maximum Length in pixels is :"<<maxLength<<endl;
	//Showing result

	imshow("Undistorted",Undistorted);
	waitKey(1);

}



void ProjectorCamera (Mat &Dst, Mat& Src, Point2f srcpt[4], Point2f dstpt[4] )
{
	Mat DGray, DBin, SGray, Window;
	vector<Vec2f> linesPC,linesRect;
	Vec2f topEdgePC, bottomEdgePC, leftEdgePC, rightEdgePC;
	Vec2f topEdgeRect, bottomEdgeRect, leftEdgeRect, rightEdgeRect;
	Point2f ptTopLeftPC, ptTopRightPC, ptBotLeftPC, ptBotRightPC;
	Point2f ptTopLeftRect, ptTopRightRect, ptBotLeftRect, ptBotRightRect;
	double alphaPC=2.2;	//1.0-3.0
	int betaPC=50;		//0-100
	Scalar mymeanPC, stddevPC;
	double stdvPC;
	
	
	if(Src.empty() || Dst.empty() )
    {
       cout<<"Image Loading failed in ProjectorCamera();"<<endl;
       getchar();
    }
	imwrite("./Original.png",Dst);
	cvtColor(Dst,DGray,CV_BGR2GRAY);
	meanStdDev(DGray, mymeanPC, stddevPC);
	stdvPC= stddevPC.val[0];
	if(stdvPC<20)
	{
		DGray.convertTo(DGray, -1, alphaPC, betaPC);
	}
	//double Imax, Imin;
	//cv::minMaxLoc(DGray, &Imin, &Imax);
	//cout<<"Imax ="<<Imax<<endl;
	blur(DGray, DGray, Size(7,7));
	threshold(DGray,DGray,235,255,THRESH_BINARY);
	//DGray=DGray>150;
	imshow("DGray",DGray);
	imwrite("./DGray.png",DGray);
	blur(DGray, DGray, Size(7,7));
	DBin=Mat(DGray.size(),CV_8UC1);
	adaptiveThreshold(DGray,DBin,200,ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,5,2);
	bitwise_not(DBin,DBin);
	dilate(DBin,DBin,Window);
	dilate(DBin,DBin,Window);
	erode(DBin,DBin,Window);
	floodFilling(DBin);
	//erode(DBin,DBin,Window);
	erode(DBin,DBin,Window);
	blur(DGray, DGray, Size(7,7));
	Canny(DBin,DBin,50,250,5);
	imshow("Dst Bin",DBin);
	imwrite("./Dst Bin.png",DBin);
	waitKey(1);


	//Rect
	SGray=Mat(Src.size(),CV_8UC1);
	cvtColor(Src,SGray,CV_BGR2GRAY);
	Canny(SGray,SGray,50,250,3);   
	//-----------------------------------------------------
	Mat NewDBin =Mat::zeros(DBin.size(),DBin.type());
	vector<vector<Point>> contours;
	vector<Vec4i> hier;
	int largestarea=0; int largestindex=0;
	findContours(DBin,contours,hier,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE
	
	for( int i = 0; i < contours.size(); i++ )
	{ 
		double iter = contourArea(contours[i],false);
		if (iter>largestarea)
		{
			largestarea=iter;
			largestindex=i;
		}	 					
	}

	vector<vector<Point>>hull( contours.size() );
	convexHull(Mat(contours[largestindex]), hull[0],false);

	Scalar color1 = Scalar(255,0,0);

	//drawContours(NewDBin, hull, 0, color1, 1, 8, vector<Vec4i>(), 0, Point());
	vector<vector<Point>>polyDp( contours.size() );

	approxPolyDP(Mat(contours[largestindex]),polyDp[0],5,true);
	//while(polyDp.size()!=4)
	//{
	//	approxPolyDP(Mat(polyDp[0]),polyDp[0],5,true);
	//}

	drawContours(NewDBin, polyDp, 0, color1, 1, 8, vector<Vec4i>(), 0, Point() );
	blur(NewDBin,NewDBin,Size(3,3));
	imshow("NewDBin",NewDBin);
	imwrite("./NewDBin.png",NewDBin);
	waitKey(1);

	//RotatedRect NewDs =minAreaRect(Mat(hull[0]));
	//Point2f NewDpts[4];
	//NewDs.points(NewDpts);
	//sortCorners(NewDpts,NewDBin.cols,NewDBin.rows);

	//ptTopLeftPC=NewDpts[0];	ptTopRightPC=NewDpts[1];
	//ptBotLeftPC=NewDpts[3];	ptBotRightPC=NewDpts[2];
	
	HoughLines(NewDBin,linesPC,1,CV_PI/180,90);
	mergeRelatedLinesPC(&linesPC,Dst);
	topEdgePC=Vec2f(1000,1000);  bottomEdgePC=Vec2f(-1000,-1000);
	leftEdgePC=Vec2f(1000,1000); rightEdgePC=Vec2f(-1000,-1000);

	HoughLines(SGray,linesRect,1,CV_PI/180,200);
	topEdgeRect=Vec2f(1000,1000); bottomEdgeRect=Vec2f(-1000,-1000);
	leftEdgeRect=Vec2f(1000,1000); rightEdgeRect=Vec2f(-1000,-1000);
    //----------------------------------------------------- 
    //Drawing PC1
	drawlines(linesPC,topEdgePC,bottomEdgePC,rightEdgePC,leftEdgePC);
    intersectionPoints(Dst,topEdgePC,bottomEdgePC,rightEdgePC,leftEdgePC,ptTopLeftPC,ptTopRightPC,ptBotLeftPC,ptBotRightPC);

	//Drawing Rect
	drawlines(linesRect,topEdgeRect,bottomEdgeRect,rightEdgeRect,leftEdgeRect);
    intersectionPoints(Src,topEdgeRect,bottomEdgeRect,rightEdgeRect,leftEdgeRect,ptTopLeftRect,ptTopRightRect,ptBotLeftRect,ptBotRightRect);
	//-----------------------------------------------------
	drawLine(topEdgePC, Dst, Scalar(0,0,0));
	drawLine(bottomEdgePC, Dst, Scalar(0,0,0));
	drawLine(rightEdgePC, Dst, Scalar(0,0,0));
	drawLine(leftEdgePC, Dst, Scalar(0,0,0));
	
	circle(Dst, ptTopLeftPC, 10, Scalar(255,0,0));
	circle(Dst, ptTopRightPC, 10, Scalar(0,255,0));
	circle(Dst, ptBotLeftPC, 10, Scalar(0,0,255));             
	circle(Dst, ptBotRightPC, 10, Scalar(255,255,0));
	imshow("Dst",Dst);
	imwrite("./Dst.png",Dst);
	
	//circle(Src, ptTopLeftRect, 10, Scalar(255,0,0));
	//circle(Src, ptTopRightRect, 10, Scalar(0,255,0));
	//circle(Src, ptBotLeftRect, 10, Scalar(0,0,255));
	//circle(Src, ptBotRightRect, 10, Scalar(255,255,0));
	//imshow("Src",Src);

	srcpt[0]=ptTopLeftRect; dstpt[0]=ptTopLeftPC;
	srcpt[1]=ptTopRightRect; dstpt[1]=ptTopRightPC;
	srcpt[2]=ptBotLeftRect; dstpt[2]=ptBotLeftPC;
	srcpt[3]=ptBotRightRect; dstpt[3]=ptBotRightPC;

	waitKey(0);// wait after every capture and conclusion
	
	Src.release();
	DGray.release();
	DBin.release();
	SGray.release();
	//destroyWindow("Dst");
	//destroyWindow("Dst Bin");
	
	//warpPerspective(TestImg,OutBlah,T,TestImg.size());
	//Showing result
	//imshow("OUTBLAH",OutBlah);
}

Mat getTransform2(Point2f srcpt[16], Point2f dstpt[16])
{
	Mat T = getPerspectiveTransform(srcpt,dstpt);
	Tinv = getPerspectiveTransform(dstpt,srcpt);
	Mat D= (Mat_<double>(16,2)<< dstpt[0].x,dstpt[0].y, dstpt[1].x,dstpt[1].y, dstpt[2].x,dstpt[2].y, dstpt[3].x,dstpt[3].y,
								dstpt[4].x,dstpt[4].y, dstpt[5].x,dstpt[5].y, dstpt[6].x,dstpt[6].y, dstpt[7].x,dstpt[7].y,
								dstpt[8].x,dstpt[8].y, dstpt[9].x,dstpt[9].y, dstpt[10].x,dstpt[10].y, dstpt[11].x,dstpt[11].y,
								dstpt[12].x,dstpt[12].y, dstpt[13].x,dstpt[13].y, dstpt[14].x,dstpt[14].y, dstpt[15].x,dstpt[15].y);
	Mat S= (Mat_<double>(16,2)<< srcpt[0].x,srcpt[0].y, srcpt[1].x,srcpt[1].y, srcpt[2].x,srcpt[2].y, srcpt[3].x,srcpt[3].y,
								srcpt[4].x,srcpt[4].y, srcpt[5].x,dstpt[5].y, srcpt[6].x,srcpt[6].y, srcpt[7].x,srcpt[7].y,
								srcpt[8].x,srcpt[8].y, srcpt[9].x,srcpt[9].y, srcpt[10].x,srcpt[10].y, srcpt[11].x,srcpt[11].y,
								srcpt[12].x,srcpt[12].y, srcpt[13].x,srcpt[13].y, srcpt[14].x,srcpt[14].y, srcpt[15].x,srcpt[15].y);
	H = findHomography(D, S, CV_RANSAC);
	cout<<"H ="<<H<<endl;
	return T;
}

Mat getTransform(Point2f srcpt[4], Point2f dstpt[4])
{
	Mat T = getPerspectiveTransform(srcpt,dstpt);
	Tinv = getPerspectiveTransform(dstpt,srcpt);
	Mat D= (Mat_<double>(4,2)<<dstpt[0].x,dstpt[0].y, dstpt[1].x,dstpt[1].y,dstpt[2].x,dstpt[2].y,dstpt[3].x,dstpt[3].y);
	Mat S= (Mat_<double>(4,2)<<srcpt[0].x,srcpt[0].y, srcpt[1].x,srcpt[1].y,srcpt[2].x,srcpt[2].y,srcpt[3].x,srcpt[3].y);
	H = findHomography(D, S, CV_RANSAC);
	cout<<"H ="<<H<<endl;
	return T;
}

Mat DetectBorder(Mat Rundistorted)
{
	Mat any, lower_bound, upper_bound, blah1;
	Mat  RGray, RoughROI, RoughROI2, RoughROI3;
	cvtColor(Rundistorted,any,CV_BGR2HSV);
	lower_bound=(Mat_<double>(1,3)<<Hu-(0.2*Hu),Sa-(0.4*Sa),Va-(0.4*Va));
	upper_bound=(Mat_<double>(1,3)<<Hu+(0.6*Hu),Sa+(0.4*Sa),Va+(0.4*Va));
	inRange(any,lower_bound,upper_bound,any);
	dilate(any,any,Window);

	imwrite("./hsv_thresholding.png",any);

	cvtColor(Rundistorted,RGray,CV_BGR2GRAY);
	//imshow("Gray",RGray);
	//waitKey(1);
	double alpha=2.2;	//1.0-3.0
	int beta=50;		//0-100
	Scalar mymean, stddev;
	double stdv;

	meanStdDev(RGray, mymean, stddev);
	stdv= stddev.val[0];
	if(stdv<20)
	{
		RGray.convertTo(RGray, -1, alpha, beta);
	}
	imwrite("./gray_equalised.png",RGray);

    GaussianBlur(RGray, RGray, Size(5,5),0);
	RoughROI = Mat(RGray.size(),CV_8UC1);
	RoughROI2 = Mat(RGray.size(),CV_8UC1);
    adaptiveThreshold(RGray,RoughROI,255,ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,5,2);
    bitwise_not(RoughROI,RoughROI);
	imwrite("./adaptive_thresholding_gray.png",RoughROI);

    Window=(Mat_<uchar>(3,3)<<0,1,0,1,1,1,0,1,0);
	
    dilate(RoughROI,RoughROI,Window);
	for (int i=0; i<10; i++)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(i)-=RoughROI.col(i);
	}
	for(int i=RoughROI.rows-1,j=RoughROI.cols-1;i>RoughROI.rows-1-20,j>RoughROI.cols-1-20;i--,j--)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(j)-=RoughROI.col(j);
	}
	imwrite("./removing_boundary_noise_gray.png",RoughROI);
	RoughROI.copyTo(blah1);
	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( any.at<uchar>(i,j)!=0){
				RoughROI.at<uchar>(i,j)=0;
			}
		}
	}
	erode(RoughROI,RoughROI,Window);
	imwrite("./NandOperation_Hsv_Gray.png",RoughROI);

	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( RGray.at<uchar>(i,j)>=180){
				blah1.at<uchar>(i,j)=0;
			}
		}
	}
	erode(blah1,blah1,Window);
	erode(blah1,blah1,Window);
		
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./thresholdingGray_extremewhite_andedwith_nandedoutput_blah1.png",blah1);
	threshold(RGray,RoughROI2,170,255,THRESH_BINARY);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI,RoughROI,Window);
	imwrite("./extremewhiteinGray_andedwithblah1_roughroi2.png",RoughROI2);
	bitwise_and(RoughROI2,RoughROI,RoughROI2);
	blah1=blah1-(RoughROI2);
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./RoughROI_andedwithroughroi2_roughROI.png",RoughROI);


	floodFilling(RoughROI);
	imwrite("./floodfilled_roughROI.png",RoughROI);

	for(int i=0;i<10;i++)
	{
		dilate(RoughROI,RoughROI,Window);
		erode(RoughROI,RoughROI,Window);
	}
	threshold(RGray,RoughROI3,20,255,THRESH_BINARY_INV);// changed 50 to 20 on experiment of Oct 5th
	for (int i=0; i<10; i++)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(i)-=RoughROI3.col(i);
	}
	for(int i=RoughROI3.rows-1,j=RoughROI3.cols-1;i>RoughROI3.rows-1-20,j>RoughROI3.cols-1-20;i--,j--)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(j)-=RoughROI3.col(j);
	}
	bitwise_or(RoughROI3,RoughROI,RoughROI);
	imwrite("./ORed_with_extreme_black_inGray.png",RoughROI);
	//imshow("Rough2",RoughROI3);
	//waitKey(1);
	return RoughROI;
}

void GetROIinScreen(Mat &ROIScreen, Mat & CROI, Mat C)
{
	if(ROIScreen.empty())
	{
		cout<<"Image Loading failed!"<<endl;
		getchar();
	}
	Mat RoughROI, Rundistorted, UndistortedROI;
	
	imwrite("./taken_frame.png", ROIScreen);
	cv::warpPerspective(ROIScreen,Rundistorted,C,Size(maxLength,floor(maxLength*0.75)));
	imwrite("./undistorted_frame.png", Rundistorted);
	

	RoughROI=DetectBorder(Rundistorted);

	floodFilling(RoughROI);
	imwrite("./Final_ROI.png",RoughROI);
	imshow("Rough",RoughROI);
	waitKey(1);
	
	GaussianBlur(RoughROI,RoughROI,Size(5,5),0);
	Canny(RoughROI,RoughROI,50,250,3);
	imwrite("./Edges_of_ROI.png",RoughROI);
	imshow("Edges",RoughROI);
	waitKey(1);
						 
	vector<vector<Point>> contours;
	vector<Vec4i> hier;
	int largestarea=0; int largestindex=0;
	findContours(RoughROI,contours,hier,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE
		
	vector<vector<Point>>hull( contours.size() );
	RotatedRect recte;
						 
	for( int i = 0; i < contours.size(); i++ )
	convexHull(Mat(contours[i]), hull[i],false);
	
	for( int i = 0; i < hull.size(); i++ )
	{ 
		double iter = contourArea(hull[i],false);
		if (iter>largestarea)
		{
			largestarea=iter;
			largestindex=i;
		}	 					
	}

	vector<vector<Point>> polydp( contours.size() );
	approxPolyDP(Mat(contours[largestindex]),polydp[0],5,true);
	convexHull(Mat(polydp[0]), hull[0],false);

	recte= minAreaRect(Mat(polydp[0]));

	RNG rng;	Mat recthull(RoughROI.size(),CV_8UC3); Mat convexhull= Mat::zeros(RoughROI.size(),CV_8UC3);
	Point2f rectpoints[4];
	cvtColor(RoughROI,recthull,CV_GRAY2BGR);				 
	Scalar color1 = Scalar(255,0,0), color2 = Scalar(0,255,0), color3 = Scalar(0,0,255); 

	drawContours( convexhull, contours, largestindex, color1, 1, 8, vector<Vec4i>(), 0, Point() );
	drawContours( convexhull, hull, 0, color2, 1, 8, vector<Vec4i>(), 0, Point() );
	drawContours( convexhull, polydp, 0, color3, 1, 8, vector<Vec4i>(), 0, Point() );
	imwrite( "./ConvexHull.png", convexhull);
	imshow("All Contours",convexhull);


	recte.points(rectpoints);
	//allocateCorners(rectpoints,Rundistorted);
	//double area_rect= rectArea(srcROI); 
	double area_rect= recte.size.area();
	cout<<"Area rect ="<<area_rect<<endl;
	
	for(int j=0; j<4;j++)
	{
		line(recthull,rectpoints[j],rectpoints[(j+1)%4],color1,1,8);
		//line(Rundistorted,rectpoints[j],rectpoints[(j+1)%4],color1,1,8);
		//circle(Rundistorted,rectpoints[j],4,color1,1,8,0);
	}
	imwrite( "./Bounding box of ConvexHull.png", recthull);
	imshow( "Bounding Box of ConvexHull", recthull);
	
	rectCentre=recte.center;	Point2f rectCentrePre=rectCentre;
	rectAngle=Point2f(recte.angle,0); 
	Size2f rectSize= recte.size;

	//int point_close= distanceCheck(srcROI);
	
	cout<<"rectAngle = "<<rectAngle.x<<endl;
	if(rectSize.width < rectSize.height)
	{
		rectAngle.x= 90+rectAngle.x;
	}

	CenfilterNoise(rectCentre);

	//cout<<"aver_W = "<<average_W<<endl;
	//cout<<"aver_H = "<<average_H<<endl;


	//if(rectAngle.x <= 5  &&  rectAngle.x >= -5 )
	//{
	//	rectAngle.x=-90;
	//}
	if(area_rect>4000 && area_rect<9000)
	{
		count_area=count_area+1;

		double width_o, height_o;
		width_o=max(rectSize.width,rectSize.height);
		height_o=min(rectSize.width,rectSize.height);

		rectSize_sumW= rectSize_sumW+width_o;
		rectSize_sumH= rectSize_sumH+height_o;

		average_W=rectSize_sumW/count_area;
		average_H=rectSize_sumH/count_area;
		cout<<"rectAngle = "<<rectAngle.x<<endl;
		AngfilterNoise(rectAngle);
		cout<<"rectAngle = "<<rectAngle.x<<endl;
	
		RotatedRect corr_rect=RotatedRect(rectCentre,Size2f((float)average_W,(float)average_H),rectAngle.x);
		Point2f corr_rectpoints[4];		corr_rect.points(corr_rectpoints);
		allocateCorners(corr_rectpoints,Rundistorted);

	
		RfilterNoise(srcROI);
		
		MoveTLpt.x=srcROI[0].x; MoveTLpt.y=srcROI[0].y;
		Mat HMT= (Mat_<double>(3,1) << MoveTLpt.x,MoveTLpt.y,1);
		
		Cinv.at<double>(0,2)=0;
		Cinv.at<double>(1,2)=0;

		HMT= Cinv*HMT;
		HMT= HMT/HMT.at<double>(2,0);

		MoveTLpt.x= HMT.at<double>(0,0);
		MoveTLpt.y= HMT.at<double>(1,0);


		double diff_cen=rectCentreDiff(rectCentre,prev_rectCentre), diff_ori =rectCentreDiff(rectCentre,Point2f(0,0));
		double diff_ang=rectAngleDiff(rectAngle,prev_rectAngle), diff_ori_a =rectAngleDiff(rectAngle,Point2f(0,0));

		cout<<"diff_cen = "<<diff_cen<<endl;
		cout<<"diff_ang = "<<diff_ang<<endl;

		if( diff_cen>2||  diff_cen==diff_ori ||diff_ang>=1 || diff_ang==diff_ori_a )//diff_cen>2||  diff_cen==diff_ori ||
		{
			CROI=getPerspectiveTransform(srcROI,dstROI);

		}

		for(int j=0; j<4;j++)
		{
			line(Rundistorted,corr_rectpoints[j],corr_rectpoints[(j+1)%4],color1,1,8);
			line(Rundistorted,srcROI[j],srcROI[(j+1)%4],color2,1,8);
			circle(Rundistorted,srcROI[j],4,color2,1,8,0);
		}
	}

	imwrite("./ConvexHull_on_undistorted_frame.png",Rundistorted);
	imshow("Rundistorted2",Rundistorted);
	waitKey(1);

	prev_rectCentre=rectCentre;
	prev_rectAngle=rectAngle;

	//cv::warpPerspective(Rundistorted,UndistortedROI,CROI,Size(maxLengthROI,floor(maxLengthROI*0.70707)));
	//imshow("UndistortedROI",UndistortedROI);
	//waitKey(1);

}

Mat MakeProjection (Mat &Screen, Mat T,Mat C, bool Method)
{
	Mat P;
    if(Screen.empty())
    {
       cout<<"Image Loading failed"<<endl;
       getchar();
    }

	if(Method==true)
		CameraScreen (Screen, C);//CS
	else
		GetROIinScreen(Screen,CROI, C);//ROI
	if(Method==true)
		P=(C)*T;
	else
	{
		//cout<<"C"<<C<<endl;
		//cout<<"CROI"<<CROI<<endl;
		//cout<<"T"<<T<<endl;
	
		P=(CROI*C)*T;

	}

	return P;
}

Mat WarpImage(Mat Image, Mat P, double ScaleX, double ScaleY, bool circular, bool triangular, Mat T)
{
	Mat warped, W;
	double WA, WB, HA, HB;
	double row, col;
	double sx, sy;
	Mat TpLt, TpRt, BtLt, BtRt;
	Mat ZTpLt, ZTpRt, ZBtLt, ZBtRt;
	Mat CMx, CMy, CP_Pt;
	double Size_Width, Size_Height;
	int w, h;

	if(Image.empty())
    {
       cout<<"Image Loading failed for WarpImage"<<endl;
       getchar();
    }
	W=(P.inv());

	row=Image.size().height;
	col=Image.size().width;

	TpLt=(Mat_<double>(3,1) << 0,0,1);		TpRt=(Mat_<double>(3,1) << 0,col,1);
	BtLt=(Mat_<double>(3,1) << row,0,1);	BtRt=(Mat_<double>(3,1) << row,col,1);
	
	ZTpLt=W*TpLt;	ZTpLt=ZTpLt/ZTpLt.at<double>(2,0);
	ZTpRt=W*TpRt;	ZTpRt=ZTpRt/ZTpRt.at<double>(2,0);
	ZBtLt=W*BtLt;	ZBtLt=ZBtLt/ZBtLt.at<double>(2,0);
	ZBtRt=W*BtRt;	ZBtRt=ZBtRt/ZBtRt.at<double>(2,0);
	
	WA=sqrt(pow((ZBtRt.at<double>(0,0)-ZBtLt.at<double>(0,0)),2)+pow((ZBtRt.at<double>(1,0)-ZBtLt.at<double>(1,0)),2));
	WB=sqrt(pow((ZTpRt.at<double>(0,0)-ZTpLt.at<double>(0,0)),2)+pow((ZTpRt.at<double>(1,0)-ZTpLt.at<double>(1,0)),2));

	HA=sqrt(pow((ZTpRt.at<double>(0,0)-ZBtRt.at<double>(0,0)),2)+pow((ZTpRt.at<double>(1,0)-ZBtRt.at<double>(1,0)),2));
	HB=sqrt(pow((ZTpLt.at<double>(0,0)-ZBtLt.at<double>(0,0)),2)+pow((ZTpLt.at<double>(1,0)-ZBtLt.at<double>(1,0)),2));
	
	Size_Width= max(WA,WB);	Size_Height= max(HA,HB);
	w= ceil(Size_Width);	h= ceil( Size_Height);
	Mat sml;
	resize(Image,sml,Size(ScaleX*Image.rows,ScaleY*Image.cols));//0.52/4, 0.20/4(previous)
	//imshow("sml",sml);

	Mat Y;
	Scalar color1 = Scalar(255,0,0);
	Scalar color2 = Scalar(0,0,255);
	//P.copyTo(Y);

	P.at<double>(0,2)=0;
	P.at<double>(1,2)=0;

	warpPerspective(sml,warped,P,Size(1920,1080),WARP_INVERSE_MAP);

	Mat G;
	Tinv.copyTo(G);
	
	proPt = CMoveTLpt+ MoveTLpt;
	CP_Pt= (Mat_<double>(3,1) << proPt.x,proPt.y,1);
	CP_Pt= G*(CP_Pt);
	CP_Pt.at<double>(0,0)= CP_Pt.at<double>(0,0)/CP_Pt.at<double>(2,0);
	CP_Pt.at<double>(1,0)= CP_Pt.at<double>(1,0)/CP_Pt.at<double>(2,0);


	//Point2f mypS[4], mypD[4];
	//mypS[0].x= 249.743; mypS[0].y= 224.183;
	//mypS[1].x= 390.809; mypS[1].y= 121.993;
	//mypS[2].x= 285.866; mypS[2].y= 197.707;
	//mypS[3].x= 187.425; mypS[3].y= 209.543;
	//mypD[0].x=540; mypD[0].y=330;
	//mypD[1].x=540; mypD[1].y=130;
	//mypD[2].x=630; mypD[2].y=230;
	//mypD[3].x=350; mypD[3].y=265;

	//Mat U = getPerspectiveTransform(mypS,mypD);

	//std::vector<Point2f> prpt(2);
	//prpt[0]= cvPoint(proPt.x,proPt.y);
	////prpt[1]= cvPoint(CMoveTLpt.x,CMoveTLpt.y);
	//std::vector<Point2f> prpt2(2);
	//std::vector<Point2f> prpt3(2);
	//perspectiveTransform(prpt,prpt2,H);
	//perspectiveTransform(prpt,prpt3,U);
	cout<<"G="<<G<<endl;                   

	cout<<"Pro_pt = ["<<proPt.x<<","<<proPt.y<<"]"<<endl;
	cout<<"CP_pt = ["<<CP_Pt.at<double>(0,0)<<","<<CP_Pt.at<double>(1,0)<<"]"<<endl;


	if(circular==true && triangular==false)
	{
		//circle(warped,Point(CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0)),5,color1,5);
		translateImg(warped,CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0));
	}
	if(triangular==true && circular==false)
	{
		//circle(warped,Point(CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0)),5,color1,5);
		translateImg(warped,CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0));
	}
	else
	{
		
		//circle(warped,Point(CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0)),5,color1,5);
		translateImg(warped,CP_Pt.at<double>(0,0),CP_Pt.at<double>(1,0));
	}
	TpLt.release();		TpRt.release();		BtLt.release();		BtRt.release();
	ZTpLt.release();	ZTpRt.release();	ZBtLt.release();	ZBtRt.release();
	return warped;
}     

Mat DetectCircle(Mat Crop)
{
	Mat Cir,Rundistorted, any, lower_bound, upper_bound, blah1;
	Mat  RGray, RoughROI, RoughROI2, RoughROI3;;
	Crop.copyTo(Rundistorted);

	cvtColor(Rundistorted,any,CV_BGR2HSV);
	lower_bound=(Mat_<double>(1,3)<<Hu-(0.2*Hu),Sa-(0.4*Sa),Va-(0.4*Va));
	upper_bound=(Mat_<double>(1,3)<<Hu+(0.6*Hu),Sa+(0.4*Sa),Va+(0.4*Va));
	inRange(any,lower_bound,upper_bound,any);
	dilate(any,any,Window);

	imwrite("./hsv_thresholding.png",any);

	cvtColor(Rundistorted,RGray,CV_BGR2GRAY);
	//imshow("Gray",RGray);
	//waitKey(1);
	double alpha=2.2;	//1.0-3.0
	int beta=50;		//0-100
	Scalar mymean, stddev;
	double stdv;

	meanStdDev(RGray, mymean, stddev);
	stdv= stddev.val[0];
	if(stdv<20)
	{
		RGray.convertTo(RGray, -1, alpha, beta);
	}
	imwrite("./gray_equalised.png",RGray);

    GaussianBlur(RGray, RGray, Size(5,5),0);
	RoughROI = Mat(RGray.size(),CV_8UC1);
	RoughROI2 = Mat(RGray.size(),CV_8UC1);
    adaptiveThreshold(RGray,RoughROI,255,ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,5,2);
    bitwise_not(RoughROI,RoughROI);
	imwrite("./adaptive_thresholding_gray.png",RoughROI);

    Window=(Mat_<uchar>(3,3)<<0,1,0,1,1,1,0,1,0);
	
    dilate(RoughROI,RoughROI,Window);
	for (int i=0; i<10; i++)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(i)-=RoughROI.col(i);
	}
	for(int i=RoughROI.rows-1,j=RoughROI.cols-1;i>RoughROI.rows-1-20,j>RoughROI.cols-1-20;i--,j--)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(j)-=RoughROI.col(j);
	}
	imwrite("./removing_boundary_noise_gray.png",RoughROI);
	RoughROI.copyTo(blah1);
	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( any.at<uchar>(i,j)!=0){
				RoughROI.at<uchar>(i,j)=0;
			}
		}
	}
	erode(RoughROI,RoughROI,Window);
	imwrite("./NandOperation_Hsv_Gray.png",RoughROI);

	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( RGray.at<uchar>(i,j)>=180){
				blah1.at<uchar>(i,j)=0;
			}
		}
	}
	erode(blah1,blah1,Window);
	erode(blah1,blah1,Window);
		
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./thresholdingGray_extremewhite_andedwith_nandedoutput_blah1.png",blah1);
	threshold(RGray,RoughROI2,170,255,THRESH_BINARY);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI,RoughROI,Window);
	imwrite("./extremewhiteinGray_andedwithblah1_roughroi2.png",RoughROI2);
	bitwise_and(RoughROI2,RoughROI,RoughROI2);
	blah1=blah1-(RoughROI2);
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./RoughROI_andedwithroughroi2_roughROI.png",RoughROI);


	floodFilling(RoughROI);
	imwrite("./floodfilled_roughROI.png",RoughROI);

	for(int i=0;i<10;i++)
	{
		dilate(RoughROI,RoughROI,Window);
		erode(RoughROI,RoughROI,Window);
	}
	threshold(RGray,RoughROI3,50,255,THRESH_BINARY_INV);
	for (int i=0; i<10; i++)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(i)-=RoughROI3.col(i);
	}
	for(int i=RoughROI3.rows-1,j=RoughROI3.cols-1;i>RoughROI3.rows-1-20,j>RoughROI3.cols-1-20;i--,j--)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(j)-=RoughROI3.col(j);
	}
	bitwise_or(RoughROI3,RoughROI,RoughROI);
	imwrite("./ORed_with_extreme_black_inGray.png",RoughROI);
	//imshow("Rough2",RoughROI3);
	//waitKey(1);

	RoughROI.copyTo(Cir);

	return Cir;
}


Mat DetectTriangle(Mat Crop)
{
	Mat Cir,Rundistorted, any, lower_bound, upper_bound, blah1;
	Mat  RGray, RoughROI, RoughROI2, RoughROI3;;
	Crop.copyTo(Rundistorted);

	cvtColor(Rundistorted,any,CV_BGR2HSV);
	lower_bound=(Mat_<double>(1,3)<<Hu-(0.2*Hu),Sa-(0.4*Sa),Va-(0.4*Va));
	upper_bound=(Mat_<double>(1,3)<<Hu+(0.6*Hu),Sa+(0.4*Sa),Va+(0.4*Va));
	inRange(any,lower_bound,upper_bound,any);
	dilate(any,any,Window);

	imwrite("./hsv_thresholding.png",any);

	cvtColor(Rundistorted,RGray,CV_BGR2GRAY);
	//imshow("Gray",RGray);
	//waitKey(1);
	double alpha=2.2;	//1.0-3.0
	int beta=50;		//0-100
	Scalar mymean, stddev;
	double stdv;

	meanStdDev(RGray, mymean, stddev);
	stdv= stddev.val[0];
	if(stdv<20)
	{
		RGray.convertTo(RGray, -1, alpha, beta);
	}
	imwrite("./gray_equalised.png",RGray);

    GaussianBlur(RGray, RGray, Size(5,5),0);
	RoughROI = Mat(RGray.size(),CV_8UC1);
	RoughROI2 = Mat(RGray.size(),CV_8UC1);
    adaptiveThreshold(RGray,RoughROI,255,ADAPTIVE_THRESH_MEAN_C,THRESH_BINARY,5,2);
    bitwise_not(RoughROI,RoughROI);
	imwrite("./adaptive_thresholding_gray.png",RoughROI);

    Window=(Mat_<uchar>(3,3)<<0,1,0,1,1,1,0,1,0);
	
    dilate(RoughROI,RoughROI,Window);
	for (int i=0; i<10; i++)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(i)-=RoughROI.col(i);
	}
	for(int i=RoughROI.rows-1,j=RoughROI.cols-1;i>RoughROI.rows-1-20,j>RoughROI.cols-1-20;i--,j--)
	{
		RoughROI.row(i)-=RoughROI.row(i);
		RoughROI.col(j)-=RoughROI.col(j);
	}
	imwrite("./removing_boundary_noise_gray.png",RoughROI);
	RoughROI.copyTo(blah1);
	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( any.at<uchar>(i,j)!=0){
				RoughROI.at<uchar>(i,j)=0;
			}
		}
	}
	erode(RoughROI,RoughROI,Window);
	imwrite("./NandOperation_Hsv_Gray.png",RoughROI);

	for(int i=0; i<RoughROI.rows-1;i++){
		for(int j=0; j<RoughROI.cols-1;j++){
			if( RGray.at<uchar>(i,j)>=180){
				blah1.at<uchar>(i,j)=0;
			}
		}
	}
	erode(blah1,blah1,Window);
	erode(blah1,blah1,Window);
		
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./thresholdingGray_extremewhite_andedwith_nandedoutput_blah1.png",blah1);
	threshold(RGray,RoughROI2,170,255,THRESH_BINARY);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI2,RoughROI2,Window);
	dilate(RoughROI,RoughROI,Window);
	imwrite("./extremewhiteinGray_andedwithblah1_roughroi2.png",RoughROI2);
	bitwise_and(RoughROI2,RoughROI,RoughROI2);
	blah1=blah1-(RoughROI2);
	bitwise_and(blah1,RoughROI,RoughROI);
	imwrite("./RoughROI_andedwithroughroi2_roughROI.png",RoughROI);


	floodFilling(RoughROI);
	imwrite("./floodfilled_roughROI.png",RoughROI);

	for(int i=0;i<10;i++)
	{
		dilate(RoughROI,RoughROI,Window);
		erode(RoughROI,RoughROI,Window);
	}
	threshold(RGray,RoughROI3,40,255,THRESH_BINARY_INV);
	for (int i=0; i<10; i++)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(i)-=RoughROI3.col(i);
	}
	for(int i=RoughROI3.rows-1,j=RoughROI3.cols-1;i>RoughROI3.rows-1-20,j>RoughROI3.cols-1-20;i--,j--)
	{
		RoughROI3.row(i)-=RoughROI3.row(i);
		RoughROI3.col(j)-=RoughROI3.col(j);
	}
	bitwise_or(RoughROI3,RoughROI,RoughROI);
	imwrite("./ORed_with_extreme_black_inGray.png",RoughROI);
	//imshow("Rough2",RoughROI3);
	//waitKey(1);

	RoughROI.copyTo(Cir);

	return Cir;
}


Mat MakeCircleProjection(Mat Screen, Mat C)
{
	Mat P, Crop, Cir, CirR;
	if(Screen.empty())
    {
       cout<<"Image Loading failed"<<endl;
       getchar();
    }
	cv::warpPerspective(Screen,Crop,C,Size(maxLength,floor(maxLength*0.75)));
	imshow("Crop",Crop);
	imwrite("./Crop_U.png",Crop);

	Cir= DetectCircle(Crop);

	floodFilling(Cir);
	imwrite("./Final_ROI_Circle.png",Cir);
	imshow("Cir",Cir);
	waitKey(1);
	
	GaussianBlur(Cir,Cir,Size(5,5),0);
	Canny(Cir,Cir,50,250,3);
	imwrite("./Edges_of_Cir.png",Cir);
	imshow("Edges",Cir);
	waitKey(1);
						 
	vector<vector<Point>> contours;
	vector<Vec4i> hier;
	int largestarea=0; int largestindex=0;
	findContours(Cir,contours,hier,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE
		
	vector<vector<Point>>hull( contours.size() );
	Rect box;
	Scalar color1 = Scalar(255,0,0), color2 = Scalar(0,255,0), color3 = Scalar(0,0,255); 

	//for( int i = 0; i < contours.size(); i++ )
	//convexHull(Mat(contours[i]), hull[i],false); not Calculating convex hull in this approach
	
	for( int i = 0; i < contours.size(); i++ )
	{ 
		double iter = contourArea(contours[i],false);
		if (iter>largestarea)
		{
			largestarea=iter;
			largestindex=i;
		}	 					
	}

	box=boundingRect(contours[largestindex]);
	cvtColor(Cir,CirR,CV_GRAY2BGR);
	drawContours(CirR,contours,largestindex,color1,1,8,vector<Vec4i>(), 0, Point() );
	rectangle( CirR,box.tl(),box.br(), color2, 2, 8, 0 );
	imshow("CirR",CirR);
	imwrite("./Circle&contours.png",CirR);


	Point2f rectPoint(box.tl()), rectPointPre=rectPoint;
	double area_rect= box.area(); 
	
	count_ave=count_ave+1;
	sum_area=sum_area+area_rect;  
	average_area=(sum_area/count_ave);
	double diff_area= abs(average_area-area_rect);
	
	sum_width=sum_width+ box.width;
	average_width=(sum_width/count_ave);

	sum_height=sum_height+ box.height;
	average_height=(sum_height/count_ave);

	cout<<"diff_area = "<< diff_area<<endl;
	cout<<"average_area = "<< average_area<<endl;
	cout<<"area_rect = "<< area_rect<<endl;
	
	bool check_area;
	if(diff_area==area_rect || diff_area<350)
		check_area=true;
	else
		check_area=false;
	if(check_area==true)
	{
		
		CenfilterNoise(rectPoint);

		Rect corr_rect=Rect(rectPoint.x,rectPoint.y,average_width,average_height);

		Point2f corr_rectpoints[4];		
		corr_rectpoints[0]=corr_rect.tl();
		corr_rectpoints[1]=corr_rect.tl()+Point(average_width,0);
		corr_rectpoints[2]=corr_rect.tl()+Point(0,average_height);
		corr_rectpoints[3]=corr_rect.br();

		allocateCorners(corr_rectpoints,Crop);
		
		MoveTLpt=srcROI[0];
		double diff_point=rectCentreDiff(rectPoint,prev_rectPoint), diff_ori_point =rectCentreDiff(rectPoint,Point2f(0,0));	
	
		CirROI=getPerspectiveTransform(srcROI,dstROI);

		for(int j=0; j<4;j++)
		{
			line(Crop,corr_rectpoints[j],corr_rectpoints[(j+1)%4],color1,1,8);
			line(Crop,srcROI[j],srcROI[(j+1)%4],color2,1,8);
			circle(Crop,srcROI[j],4,color2,1,8,0);
		}
	}

	imwrite("./ConvexHull_on_Crop.png",Crop);
	imshow("Full_Crop",Crop);
	waitKey(1);
    prev_area=area_rect;
	prev_rectCentre=rectCentre;
	prev_rectAngle=rectAngle;

	return CirROI;
}

Mat Make_n_Warp(Mat Screen, Mat T, Mat C, Mat pic)
{
	Mat Pc, Output;
	Pc= MakeCircleProjection(Screen, C);
	Output=WarpImage(pic,Pc,0.17,0.12, true,false, T);

	return Output;

}

Mat MakeTriProjection(Mat Screen, Mat C)
{
	Mat Cropped, Tri, TriR;
	if(Screen.empty())
    {
       cout<<"Image Loading failed"<<endl;
       getchar();
    }
	cv::warpPerspective(Screen,Cropped,C,Size(maxLength,floor(maxLength*0.75)));
	imshow("Crop",Cropped);
	imwrite("./Crop_U.png",Cropped);
	waitKey(1);

	Tri= DetectTriangle(Cropped);
	imshow("Tri", Tri);
	waitKey(1);	

	floodFilling(Tri);
	imwrite("./Final_ROI_Tri.png",Tri);
	imshow("Tri",Tri);
	waitKey(1);
	
	GaussianBlur(Tri,Tri,Size(5,5),0);
	Canny(Tri,Tri,50,250,3);
	imwrite("./Edges_of_Tri.png",Tri);
	imshow("Edges",Tri);
	waitKey(1);

	vector<vector<Point>> contours;
	vector<Vec4i> hier;
	int largestarea=0; int largestindex=0;
	findContours(Tri,contours,hier,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE
		
	vector<vector<Point>>hull( contours.size() );
	Rect box;
	Scalar color1 = Scalar(255,0,0), color2 = Scalar(0,255,0), color3 = Scalar(0,0,255); 

	//for( int i = 0; i < contours.size(); i++ )
	//convexHull(Mat(contours[i]), hull[i],false); not Calculating convex hull in this approach
	
	for( int i = 0; i < contours.size(); i++ )
	{ 
		double iter = contourArea(contours[i],false);
		if (iter>largestarea)
		{
			largestarea=iter;
			largestindex=i;
		}	 					
	}

	box= boundingRect(contours[largestindex]);
	cvtColor(Tri,TriR,CV_GRAY2BGR);
	drawContours(TriR,contours,largestindex,color1,1,8,vector<Vec4i>(), 0, Point() );
	rectangle( TriR,box.tl(),box.br(), color2, 2, 8, 0 );
	imshow("TriR",TriR);
	imwrite("./Triangle&contours.png",TriR);


	Point2f rectPoint(box.tl()), rectPointPre=rectPoint;
	double area_rect= box.area(); 
	
	count_ave=count_ave+1;
	sum_area=sum_area+area_rect;  
	average_area=(sum_area/count_ave);
	double diff_area= abs(average_area-area_rect);
	
	sum_width=sum_width+ box.width;
	average_width=(sum_width/count_ave);

	sum_height=sum_height+ box.height;
	average_height=(sum_height/count_ave);

	cout<<"diff_area = "<< diff_area<<endl;
	cout<<"average_area = "<< average_area<<endl;
	cout<<"area_rect = "<< area_rect<<endl;
	
	bool check_area;
	if(diff_area==area_rect || diff_area<350)
		check_area=true;
	else
		check_area=false;
	if(check_area==true)
	{
		
		CenfilterNoise(rectPoint);

		Rect corr_rect=Rect(rectPoint.x,rectPoint.y,average_width,average_height);

		Point2f corr_rectpoints[4];		
		corr_rectpoints[0]=corr_rect.tl();
		corr_rectpoints[1]=corr_rect.tl()+Point(average_width,0);
		corr_rectpoints[2]=corr_rect.tl()+Point(0,average_height);
		corr_rectpoints[3]=corr_rect.br();

		allocateCorners(corr_rectpoints,Cropped);
		
		MoveTLpt=srcROI[0];
		double diff_point=rectCentreDiff(rectPoint,prev_rectPoint), diff_ori_point =rectCentreDiff(rectPoint,Point2f(0,0));	
	
		CTriROI=getPerspectiveTransform(srcROI,dstROI);

		for(int j=0; j<4;j++)
		{
			line(Cropped,corr_rectpoints[j],corr_rectpoints[(j+1)%4],color1,1,8);
			line(Cropped,srcROI[j],srcROI[(j+1)%4],color2,1,8);
			circle(Cropped,srcROI[j],4,color2,1,8,0);
		}
	}

	imwrite("./ConvexHull_on_Crop.png",Cropped);
	imshow("Full_Crop",Cropped);
	waitKey(1);
    prev_area=area_rect;
	prev_rectCentre=rectCentre;
	prev_rectAngle=rectAngle;

	return CTriROI;
}

Mat Tri_Make_n_Warp(Mat Screen, Mat T, Mat C, Mat pic)
{
	Mat Pc, Output;
	Pc= MakeTriProjection(Screen, C);
	Output=WarpImage(pic,Pc,0.17,0.12, false, true, T);

	return Output;

}

void SetExposure(VideoCapture Cap)
{
	Mat Fr;
	double i= Cap.get(CV_CAP_PROP_EXPOSURE);
	cout<<"i= "<<i<<endl;
	while(1)
	{
		Cap>>Fr;
		imshow("SET EXPOSURE",Fr);
		waitKey(1);
		if(GetAsyncKeyState(VK_UP))
		{
			if(i<=-20)
				i=-20;
			else
				i+=0.5;
		}
		if(GetAsyncKeyState(VK_DOWN))
		{
			if(i>=20)
				i=20;
			else
				i-=0.5;
		}
		Cap.set(CV_CAP_PROP_EXPOSURE,i);
		if (GetAsyncKeyState(VK_SPACE))
			break;
	}
	Fr.release();
	destroyWindow("SET EXPOSURE");

}

void ShowSrcImage(Mat Src)
{
	GetDesktopResolution(ho,ve); 
	namedWindow("Src Image",1);
	moveWindow("Src Image",(ho-15),(0-30));//(ho-15),(0-30)
	//resize(Src,Src,Size(1920,1080));
	//resizeWindow("Projected", wid-5, hei-5);
	imshow("Src Image",Src);
	waitKey(1);
}

Mat Capturing(Mat Src)
{
	ShowSrcImage(Src);
	Mat Required, Ro;
	again:
	VideoCapture Cap(0);
	if(!Cap.isOpened())
	{	
		cout<<" Can't open camera for projector"<<endl;
		goto again;
	}


	VideoWriter SrcPtFinding;
	Size S = Size((int) Cap.get(CV_CAP_PROP_FRAME_WIDTH),(int) Cap.get(CV_CAP_PROP_FRAME_HEIGHT));
	SrcPtFinding.open("SrcPtFinding.avi",CV_FOURCC('D','I','V','4'), 60, S, true);
	if(!SrcPtFinding.isOpened())
		cout<<"Can't save video"<<endl;
	
	SetExposure(Cap);
	cout<<"Press ENTER to capture picture from camera output"<<endl;
	while(1)
	{
		Cap>>Ro;
		SrcPtFinding<<Ro;
		imshow("Press Enter to Select",Ro);
		if(waitKey(10)==13)
			break;
	}
	Ro.copyTo(Required);
	imshow("Selected Frame",Required);
	Cap.release();
	Ro.release();
	destroyWindow("Press Enter to Select");
	destroyWindow("Projected");
	destroyWindow("Src Image");
	return Required;
}

void on_H(int, void*)
{
	Hu=(double)H_slide;
}

void on_S(int, void*)
{
	Sa=(double)S_slide;
}

void on_V(int, void*)
{
	Va=(double)V_slide;
}

Mat HSV(VideoCapture & Cap, Mat& Input)
{
	Mat Output;
	Mat lower_bound,upper_bound;
	//72.1, 100, 50.4
	if(call==true)
	{
		namedWindow("Out",1);
		createTrackbar("H","Out",&H_slide,maxH,on_H);
		createTrackbar("S","Out",&S_slide,maxS,on_S);
		createTrackbar("V","Out",&V_slide,maxV,on_V);
		cout<<"Press Space once done"<<endl;
		while(!GetAsyncKeyState(VK_SPACE))	
		{
			Cap>>Input;
			cvtColor(Input,Output,CV_BGR2HSV);
			lower_bound=(Mat_<double>(1,3)<<Hu-(0.4*Hu),Sa-(0.4*Sa),Va-(0.4*Va));
			upper_bound=(Mat_<double>(1,3)<<Hu+(0.4*Hu),Sa+(0.4*Sa),Va+(0.4*Va));
			inRange(Output,lower_bound,upper_bound,Output);
			imshow("Out",Output);
			imshow("Input",Input);
			waitKey(1);
		}
		destroyWindow("Input");
		destroyWindow("Output");
		call= false;
	}
	else
	{
		Cap>>Input;
		cvtColor(Input,Output,CV_BGR2HSV);
		lower_bound=(Mat_<double>(1,3)<<Hu-(0.4*Hu),Sa-(0.4*Sa),Va-(0.4*Va));
		upper_bound=(Mat_<double>(1,3)<<Hu+(0.4*Hu),Sa+(0.4*Sa),Va+(0.4*Va));
		inRange(Output,lower_bound,upper_bound,Output);
		imshow("Out",Output);
	}
	return Output;
}

void filterInit(Point2f s[4])
{
	KF1.transitionMatrix = *(Mat_<float>(4,4)<<1,0,1,0, 0,1,0,1, 0,0,1,0, 0,0,0,1);
	measure1.setTo(Scalar(0));
	KF1.statePre.at<float>(0) = s[0].x; KF1.statePre.at<float>(1) = s[0].y;
	KF1.statePre.at<float>(2) = 0;		KF1.statePre.at<float>(3) = 0;

	KF2.transitionMatrix = *(Mat_<float>(4,4)<<1,0,1,0, 0,1,0,1, 0,0,1,0, 0,0,0,1);
	measure2.setTo(Scalar(0));
	KF2.statePre.at<float>(0) = s[1].x;	KF2.statePre.at<float>(1) = s[1].y;
	KF2.statePre.at<float>(2) = 0;		KF2.statePre.at<float>(3) = 0;

	KF3.transitionMatrix = *(Mat_<float>(4,4)<<1,0,1,0, 0,1,0,1, 0,0,1,0, 0,0,0,1);
	measure3.setTo(Scalar(0));
	KF3.statePre.at<float>(0) = s[2].x;	KF3.statePre.at<float>(1) = s[2].y;
	KF3.statePre.at<float>(2) = 0;		KF3.statePre.at<float>(3) = 0;

	KF4.transitionMatrix = *(Mat_<float>(4,4)<<1,0,1,0, 0,1,0,1, 0,0,1,0, 0,0,0,1);
	measure4.setTo(Scalar(0));
	KF4.statePre.at<float>(0) = s[3].x;	KF4.statePre.at<float>(1) = s[3].y;
	KF4.statePre.at<float>(2) = 0;		KF4.statePre.at<float>(3) = 0;

	setIdentity(KF1.measurementMatrix);						setIdentity(KF1.processNoiseCov,Scalar::all(1e-4));
	setIdentity(KF1.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(KF1.errorCovPost,Scalar::all(.1));

	setIdentity(KF2.measurementMatrix);						setIdentity(KF2.processNoiseCov,Scalar::all(1e-4));
	setIdentity(KF2.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(KF2.errorCovPost,Scalar::all(.1));

	setIdentity(KF3.measurementMatrix);						setIdentity(KF3.processNoiseCov,Scalar::all(1e-4));
	setIdentity(KF3.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(KF3.errorCovPost,Scalar::all(.1));

	setIdentity(KF4.measurementMatrix);						setIdentity(KF4.processNoiseCov,Scalar::all(1e-4));
	setIdentity(KF4.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(KF4.errorCovPost,Scalar::all(.1));
}

void filterNoise(Point2f s[4])
{
	Mat prediction1 = KF1.predict(); Mat prediction2 = KF2.predict();
	Mat prediction3 = KF3.predict(); Mat prediction4 = KF4.predict();

	Point2f predictPt1(prediction1.at<float>(0),prediction1.at<float>(1));
	Point2f predictPt2(prediction2.at<float>(0),prediction2.at<float>(1));
	Point2f predictPt3(prediction3.at<float>(0),prediction3.at<float>(1));
	Point2f predictPt4(prediction4.at<float>(0),prediction4.at<float>(1));

	measure1(0) = s[0].x; measure1(1) = s[0].y; measure2(0) = s[1].x; measure2(1) = s[1].y;
	measure3(0) = s[2].x; measure3(1) = s[2].y; measure4(0) = s[3].x; measure4(1) = s[3].y;

	Point2f measPt1(measure1(0),measure1(1)); Point2f measPt2(measure2(0),measure2(1));
	Point2f measPt3(measure3(0),measure3(1)); Point2f measPt4(measure4(0),measure4(1));

	Mat estimated1 = KF1.correct(measure1); Mat estimated2 = KF2.correct(measure2);
	Mat estimated3 = KF3.correct(measure3); Mat estimated4 = KF4.correct(measure4);
	
	Point2f statePt1(estimated1.at<float>(0),estimated1.at<float>(1)); Point2f statePt2(estimated2.at<float>(0),estimated2.at<float>(1));
	Point2f statePt3(estimated3.at<float>(0),estimated3.at<float>(1)); Point2f statePt4(estimated4.at<float>(0),estimated4.at<float>(1));

	s[0].x = statePt1.x; s[0].y = statePt1.y;	s[1].x = statePt2.x; s[1].y = statePt2.y;
	s[2].x = statePt3.x; s[2].y = statePt3.y;	s[3].x = statePt4.x; s[3].y = statePt4.y;
}

void CenfilterInit(Point2f &s)
{
	CenKF.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Cenmeasure.setTo(Scalar(0));
	CenKF.statePre.at<float>(0) = s.x;		CenKF.statePre.at<float>(1) = s.y;
	CenKF.statePre.at<float>(2) = 0;		CenKF.statePre.at<float>(3) = 0;
	CenKF.statePre.at<float>(4) = 0;		CenKF.statePre.at<float>(5) = 0;
	CenKF.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(CenKF.processNoiseCov,Scalar::all(1e-4));
	setIdentity(CenKF.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(CenKF.errorCovPost,Scalar::all(.1));
	
}

void CenfilterNoise(Point2f &s)
{
	Mat Cenprediction = CenKF.predict(); 
	Point2f CenpredictPt(Cenprediction.at<float>(0),Cenprediction.at<float>(1));
	Cenmeasure(0) = s.x; Cenmeasure(1) = s.y;
	Point2f CenmeasPt(Cenmeasure(0),Cenmeasure(1));
	Mat Cestimated = CenKF.correct(Cenmeasure);
	Point2f CenstatePt(Cestimated.at<float>(0),Cestimated.at<float>(1));
	s.x = CenstatePt.x; s.y = CenstatePt.y;
}

void AngfilterInit(Point2f &s)
{
	AngKF.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Angmeasure.setTo(Scalar(0));
	AngKF.statePre.at<float>(0) = s.x;		AngKF.statePre.at<float>(1) = s.y;
	AngKF.statePre.at<float>(2) = 0;		AngKF.statePre.at<float>(3) = 0;
	AngKF.statePre.at<float>(4) = 0;		AngKF.statePre.at<float>(5) = 0;
	AngKF.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(AngKF.processNoiseCov,Scalar::all(1e-4));
	setIdentity(AngKF.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(AngKF.errorCovPost,Scalar::all(.1));
	
}

void AngfilterNoise(Point2f &s)
{
	Mat Angprediction = AngKF.predict(); 
	Point2f AngpredictPt(Angprediction.at<float>(0),Angprediction.at<float>(1));
	Angmeasure(0) = s.x; Angmeasure(1) = s.y;
	Point2f AngmeasPt(Angmeasure(0),Angmeasure(1));
	Mat Aestimated = AngKF.correct(Angmeasure);
	Point2f AngstatePt(Aestimated.at<float>(0),Aestimated.at<float>(1));
	s.x = AngstatePt.x; s.y = AngstatePt.y;
}

void RfilterInit(Point2f s[4])
{
	RKF1.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Rmeasure1.setTo(Scalar(0));
	RKF1.statePre.at<float>(0) = s[0].x; RKF1.statePre.at<float>(1) = s[0].y;
	RKF1.statePre.at<float>(2) = 0;		RKF1.statePre.at<float>(3) = 0;
	RKF1.statePre.at<float>(4) = 0;		RKF1.statePre.at<float>(5) = 0;

	RKF2.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Rmeasure2.setTo(Scalar(0));
	RKF2.statePre.at<float>(0) = s[1].x;	RKF2.statePre.at<float>(1) = s[1].y;
	RKF2.statePre.at<float>(2) = 0;		RKF2.statePre.at<float>(3) = 0;
	RKF2.statePre.at<float>(4) = 0;		RKF2.statePre.at<float>(5) = 0;

	RKF3.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Rmeasure3.setTo(Scalar(0));
	RKF3.statePre.at<float>(0) = s[2].x;	RKF3.statePre.at<float>(1) = s[2].y;
	RKF3.statePre.at<float>(2) = 0;		RKF3.statePre.at<float>(3) = 0;
	RKF3.statePre.at<float>(4) = 0;		RKF3.statePre.at<float>(5) = 0;

	RKF4.transitionMatrix = *(Mat_<float>(6,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5, 0,0,1,0,1,0, 0,0,0,1,0,1, 0,0,0,0,1,0, 0,0,0,0,0,1 );
	Rmeasure4.setTo(Scalar(0));
	RKF4.statePre.at<float>(0) = s[3].x;	RKF4.statePre.at<float>(1) = s[3].y;
	RKF4.statePre.at<float>(2) = 0;		RKF4.statePre.at<float>(3) = 0;
	RKF4.statePre.at<float>(4) = 0;		RKF4.statePre.at<float>(5) = 0;

	RKF1.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(RKF1.processNoiseCov,Scalar::all(1e-4));
	setIdentity(RKF1.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(RKF1.errorCovPost,Scalar::all(.1));

	RKF2.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(RKF2.processNoiseCov,Scalar::all(1e-4));
	setIdentity(RKF2.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(RKF2.errorCovPost,Scalar::all(.1));

	RKF3.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(RKF3.processNoiseCov,Scalar::all(1e-4));
	setIdentity(RKF3.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(RKF3.errorCovPost,Scalar::all(.1));

	RKF4.measurementMatrix=*(Mat_<float>(2,6)<<1,0,1,0,0.5,0, 0,1,0,1,0,0.5);
	setIdentity(RKF4.processNoiseCov,Scalar::all(1e-4));
	setIdentity(RKF4.measurementNoiseCov,Scalar::all(1e-1));	setIdentity(RKF4.errorCovPost,Scalar::all(.1));
}

void RfilterNoise(Point2f s[4])
{
	Mat Rprediction1 = RKF1.predict(); Mat Rprediction2 = RKF2.predict();	Mat Rprediction3 = RKF3.predict(); Mat Rprediction4 = RKF4.predict();

	Point2f RpredictPt1(Rprediction1.at<float>(0),Rprediction1.at<float>(1));
	Point2f RpredictPt2(Rprediction2.at<float>(0),Rprediction2.at<float>(1));
	Point2f RpredictPt3(Rprediction3.at<float>(0),Rprediction3.at<float>(1));
	Point2f RpredictPt4(Rprediction4.at<float>(0),Rprediction4.at<float>(1));

	Rmeasure1(0) = s[0].x; Rmeasure1(1) = s[0].y; Rmeasure2(0) = s[1].x; Rmeasure2(1) = s[1].y;
	Rmeasure3(0) = s[2].x; Rmeasure3(1) = s[2].y; Rmeasure4(0) = s[3].x; Rmeasure4(1) = s[3].y;

	Point2f RmeasPt1(Rmeasure1(0),Rmeasure1(1)); Point2f RmeasPt2(Rmeasure2(0),Rmeasure2(1));
	Point2f RmeasPt3(Rmeasure3(0),Rmeasure3(1)); Point2f RmeasPt4(Rmeasure4(0),Rmeasure4(1));

	Mat Restimated1 = RKF1.correct(Rmeasure1); Mat Restimated2 = RKF2.correct(Rmeasure2);
	Mat Restimated3 = RKF3.correct(Rmeasure3); Mat Restimated4 = RKF4.correct(Rmeasure4);
	
	Point2f RstatePt1(Restimated1.at<float>(0),Restimated1.at<float>(1)); Point2f RstatePt2(Restimated2.at<float>(0),Restimated2.at<float>(1));
	Point2f RstatePt3(Restimated3.at<float>(0),Restimated3.at<float>(1)); Point2f RstatePt4(Restimated4.at<float>(0),Restimated4.at<float>(1));

	s[0].x = RstatePt1.x; s[0].y = RstatePt1.y;	s[1].x = RstatePt2.x; s[1].y = RstatePt2.y;
	s[2].x = RstatePt3.x; s[2].y = RstatePt3.y;	s[3].x = RstatePt4.x; s[3].y = RstatePt4.y;
}

Mat ProCamMatrix(bool debug)
{
	Mat T;
	Mat Src, Src2, Src3, Src4;
	Mat Dst, Dst2, Dst3, Dst4;
	Point2f srcpt[4], dstpt[4];
	//Point2f srcpt2[4], dstpt2[4];
	//Point2f srcpt3[4], dstpt3[4];
	//Point2f srcpt4[4], dstpt4[4];
	//Point2f srcp[16], dstp[16];
	
	//Src = imread("black2.png",1);
	//Src2 = imread("black3.png",1);
	//Src3 = imread("black4.png",1);
	//Src4 = imread("black5.png",1);

	Src = imread("calib.png",1);
	//Src2 = imread("calib1.png",1);
	//Src3 = imread("calib2.png",1);
	//Src4 = imread("calib3.png",1);
	if(debug==true)
	{
	/*	VideoCapture  Drt;
		Drt.open("DstDebug.avi");*/
		//if(!Drt.isOpened())
		//{
		//	cout<<"DstDebug not opened!"<<endl;
		//	getchar();
		//}
		/*Drt>>Dst;*/
		//if(Dst.empty())
		//{
		//	cout<<"Dst is empty"<<endl;
		//	getchar();
		//}
		Dst = imread("DstDebug.png");
	}
	else
	{
		Dst= Capturing(Src);
		//Dst2= Capturing(Src2);
		//Dst3= Capturing(Src3);
		//Dst4= Capturing(Src4);
	}
	if(debug==true)
	{
		ProjectorCamera(Dst, Src, srcpt, dstpt); //PC
		T= getTransform(srcpt, dstpt);
	}

	else
	{
		ProjectorCamera(Dst, Src, srcpt, dstpt); //PC
		//ProjectorCamera(Dst2, Src2, srcpt2, dstpt2); //PC
		//ProjectorCamera(Dst3, Src3, srcpt3, dstpt3); //PC
		//ProjectorCamera(Dst4, Src4, srcpt4, dstpt4); //PC
		/*for (int i =0; i<4; i++)
		{
			srcp[i]=srcpt[i]; dstp[i]=dstpt[i];
			srcp[i+4]=srcpt2[i]; dstp[i+4]=dstpt2[i];
			srcp[i+8]=srcpt3[i]; dstp[i+8]=dstpt3[i];
			srcp[i+12]=srcpt4[i]; dstp[i+12]=dstpt4[i];
		}*/
		T= getTransform(srcpt,dstpt);
	}
	return T;
}

void identifyShape(Mat Screen, bool& circular, bool& triangular)
{
	if(Screen.empty())
	{
		cout<<"Image Loading failed! stopped at identifying shape"<<endl;
		getchar();
	}
	Mat RoSc, RoScrn, Cntr;
	cv::warpPerspective(Screen,RoSc,C,Size(maxLength,floor(maxLength*0.75)));	

	RoScrn=DetectBorder(RoSc);

	floodFilling(RoScrn);
	//imshow("RoScrn",RoScrn);
	//waitKey(1);
	
	GaussianBlur(RoScrn,RoScrn,Size(5,5),0);
	Canny(RoScrn,RoScrn,50,250,3);
	//imshow("Edges",RoScrn);
	//waitKey(1);
						 
	vector<vector<Point>> contours;
	vector<Vec4i> hier;
	int largestarea=0; int largestindex=0;
	findContours(RoScrn,contours,hier,CV_RETR_EXTERNAL,CV_CHAIN_APPROX_NONE,Point(0, 0));//CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE

	vector<vector<Point>>hull(contours.size());				 
	for( int i = 0; i < contours.size(); i++ )
	{ 
		double iter = contourArea(contours[i],false);
		if (iter>largestarea)
		{
			largestarea=iter;
			largestindex=i;
		}	 					
	}

	vector<vector<Point>> polydp(contours.size());
	approxPolyDP(Mat(contours[largestindex]),polydp[0],5,true);
	convexHull(Mat(polydp[0]), hull[0],false);

	Scalar color1 = Scalar(255,0,0), color2 = Scalar(0,255,0), color3 = Scalar(0,0,255); 
	Cntr=Mat::zeros(RoScrn.size(),RoScrn.type());
	drawContours( Cntr, polydp, 0, color1, 1, 8, vector<Vec4i>(), 0, Point());
	imshow("Cntr",Cntr);
	waitKey(1);

	cout<<"Hull Corners ="<<hull[0].size()<<endl;
	if(hull[0].size()==4 || hull[0].size()==5)
	{
		circular = false;
		triangular = false;
		cout<<"Rectangle detected"<<endl;
	}
	else if(hull[0].size()==3)
	{
		circular = false;
		triangular = true;
		cout<<"Triangle detected"<<endl;
	}
	else if(hull[0].size()>6)
	{
		circular = true;
		triangular = false;
		cout<<"Circle detected"<<endl;
	}
	else
	{
		cout<<"shape unsure"<<endl;
	}                

}


int main ()
{	//-----------Filter Initialization-----------//

	filterInit(src);
	CenfilterInit(rectCentre);
	AngfilterInit(rectAngle);
	RfilterInit(srcROI);
	GetDesktopResolution(ho,ve); 

	//---------Projector Camera Part-------------//
	bool debug=false;
	Mat T = ProCamMatrix(debug);


	//-----------Screen Video------------------------//
	
	circular = false;//====================>>>>> circle or not
	triangular = false;//==================>>>>> triangle or not

	VideoCapture CameraCap(0);
	
	//VideoCapture CameraCap;
	//if(circular==false && triangular== false)
	//{
	//	CameraCap.open("trial.avi");
	//	//CameraCap.open("Circle2.avi");
	//	//CameraCap.open("Tri2.avi");
	//	//CameraCap.open("22.avi");
	//	//CameraCap.open("Rect2.avi");
	//	//CameraCap.open("Rect0.avi");
	//	//CameraCap.open("This.avi");//screen.avi   ScreenVideo.mp4
	//}
	//else if (circular==true && triangular == false)
	//{
	//	CameraCap.open("Circle2.avi");
	//}
	//else
	//{
	//	CameraCap.open("Tri2.avi");
	//}

	if(!CameraCap.isOpened())
	{	
		cout<<" Can't open video of screen"<<endl;
		return -1;
	}
	//-----------------Video Writer------------------//
	VideoWriter CameraFeed;
	Size S = Size((int) CameraCap.get(CV_CAP_PROP_FRAME_WIDTH),(int) CameraCap.get(CV_CAP_PROP_FRAME_HEIGHT));
	CameraFeed.open("CameraFeed.avi",CV_FOURCC('D','I','V','4'), 60, S, true);
	if(!CameraFeed.isOpened())
		cout<<"Can't save video"<<endl;

	//Size S1 = Size((int) CameraCap.get(CV_CAP_PROP_FRAME_WIDTH),(int) CameraCap.get(CV_CAP_PROP_FRAME_HEIGHT));
	//CameraFeed1.open("CameraFeed1.avi",CV_FOURCC('D','I','V','4'), 30, S1, true);
	//if(!CameraFeed1.isOpened())
	//	cout<<"Can't save video"<<endl;
	//
	//Size S2 = Size((int) CameraCap.get(CV_CAP_PROP_FRAME_WIDTH),(int) CameraCap.get(CV_CAP_PROP_FRAME_HEIGHT));
	//CameraFeed2.open("CameraFeed2.avi",CV_FOURCC('D','I','V','4'), 30, S2, true);
	//if(!CameraFeed2.isOpened())
	//	cout<<"Can't save video"<<endl;

	//-----------------Projector Video Or Image------//
	//bool projVideo=false;
	//VideoCapture Projecting;
	//if(projVideo==true)
	//{
	//	Projecting.open("MyVideo.mp4");
	//	if(!Projecting.isOpened())
	//	{	
	//		cout<<" Can't open video"<<endl;
	//		return -1;
	//	}
	//}

	bool end=false;
	//-----------------Main Loop----------------------//
	
	Mat TakenFrame;
	Mat P, Output, Screen;

	call = true;
	method= false;//==========> change to false for ROI
	
	while(end==false)
	{
get_another:
get_next_frame:

		waitKey(1);
		//Taking Picture to Detect Screen
		CameraCap>>Screen;
		if(Screen.empty())
			break;
		//VideoWriting
		CameraFeed<<Screen;
		//Making Matrix to Warp
		if(call==true && method ==false)
		{
			CameraScreen(Screen,C);
			if(GetAsyncKeyState(VK_SHIFT))
				return -1;
			if(!GetAsyncKeyState(VK_BACK))
				goto get_next_frame;

		}
		destroyWindow("Undistorted"); destroyWindow("Border");
		call= false;

		//identifyShape(Screen, circular, triangular);//==================================> shape identification
		
		if(circular == false && triangular == false)
		{
			P = MakeProjection(Screen,T,C,method);
						
			cout<<"max length is : "<<maxLength<<endl;//========> max length couts
			if(maxLength<15 || maxLength>Screen.cols )
			{
				cout<<"maxLength out of range, skiping the loop"<<endl;
				continue;
			}
			cout<<"max length ROI is : "<<maxLengthROI<<endl;//========> max length ROI couts
			if(maxLengthROI<5 || maxLengthROI>maxLength )
			{
				cout<<"maxLengthROI out of range, skiping the loop"<<endl;
				continue;
			}
		
			TakenFrame=imread("sunset.jpg",1);
			if(TakenFrame.empty())
				end=true;
			//Warping
			//pyrDown(TakenFrame,TakenFrame,Size(TakenFrame.cols/2, TakenFrame.rows/2));
			//Output= WarpImage(TakenFrame,P,0.75,1.05);//0.6,1.0);//0.65,1.15//0.45, 0.7=======>>>>> for CameraScreen only
			//Output= WarpImage(TakenFrame,P,0.185,0.23, false);//0.1,0.3);//0.35,0.5);//================>>>>>>> for ROI 
			Output= WarpImage(TakenFrame,P,0.135,0.05, false,false, T);//0.1,0.3);//0.35,0.5);//================>>>>>>> for ROI 
			//pyrUp(Output,Output, Size(Output.cols*2, Output.rows*2 ));
			namedWindow("Displaying",1);
			imshow("Displaying",Output);  
			moveWindow("Displaying",(ho-15),(0-30));
			waitKey(1);

			circle(Screen,proPt,3,Scalar(0,255,0));
			circle(Screen,CMoveTLpt,3,Scalar(0,0,255));
			imshow("Screen",Screen);
			waitKey(1);

			CameraFeed<<Screen;

		}
		else if (triangular == false && circular==true)
		{
			Mat pic = imread("sunsetCircular.jpg",1);
			if(pic.empty())
				end=true;
			
			Output=Make_n_Warp(Screen, T, C, pic);
			//resize(Output,Output,Size(1280,800));
			imshow("Displaying",Output);
			waitKey(1);

		}
		else
		{
			
			Mat pic = imread("triangle.jpg",1); 
			if(pic.empty())
				end=true;
			
			Output=Tri_Make_n_Warp(Screen, T, C, pic);
			//resize(Output,Output,Size(1280,800));
			imshow("Displaying",Output);
			waitKey(1);
		}
		//Escape          
		if(GetAsyncKeyState(VK_SHIFT))
			break;

		//while(waitKey(0)==NULL);//================= pause per frame
		//waitKey(1);
	}
	//----------------------End---------------//
	return 0;
}
