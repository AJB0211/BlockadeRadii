#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <chrono>
#include <cstdlib> //for rand function
using namespace std;

double pi = 3.141592;
double sigma = 10;



double dist(double x1,double y1, double x2, double y2){
	double d = pow(x1-x2,2) + pow(y1-y2,2);
	return sqrt(d);
}

double randDouble(double low, double high){	
	double temp;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution (low, high);
	temp = distribution(generator);
	return temp;
}

double gaussian(double x, double sig){
	double gexp = -1*pow(x,2);
	gexp /= 2*pow(sig,2);
	double g = exp(gexp);
	//g = pow(sig*sqrt(2*pi),-2);
	return g;
}


bool distcheck(double x1, double y1, double x2, double y2, double r){
	if(dist(x1,y1,x2,y2)<r){
		//if(randDouble(0,100)<5){cout << "returned false" << endl;}
		return false;
	}
	else{
		//if(randDouble(0,100)<5){cout << "returned true" << endl;}
		return true;
	}
}

//Trapezoidal integration of a Gaussian function in 2D
//steps is half the number of total steps used
double gausstrap(int steps, double range){
	double h = range/(2*steps); //stepsize
	double xiyi = 0;
	double xic = 0;
	double xid = 0;
	double ayi = 0;
	double byi = 0;

	double corners = gaussian(h*steps,sigma)*gaussian(h*(-1)*steps,sigma);
	corners += gaussian(h*(-1)*steps,sigma)*gaussian(steps*h,sigma);
	corners += gaussian(h*(-1)*steps,sigma)*gaussian(h*(-1)*steps,sigma);
	corners += gaussian(h*steps,sigma)*gaussian(h*steps,sigma);
	
	for(int k=(-1)*steps+1;k<steps;k++){
		xic += 2*gaussian(k*h,sigma)*gaussian(-1*steps*h,sigma);
		xid += 2*gaussian(k*h,sigma)*gaussian(steps*h,sigma);
		ayi += 2*gaussian((-1)*steps*h,sigma)*gaussian(k*h,sigma);
		byi += 2*gaussian(steps*h,sigma)*gaussian(k*h,sigma);
	}


	for(int kx=-1*steps+1;kx<steps;kx++){
		for(int ky=-1*steps+1;ky<steps;ky++)
			xiyi += 4*gaussian(kx*h,sigma)*gaussian(ky*h,sigma);
	}

	//cout << "h= " << h << endl;
	return pow(h,2)*(xiyi+xic+xid+ayi+byi+corners)/(pow(sigma*sqrt(2*pi),2)*4);
}





int main(int argc, char * const argv[]){
	int imax = 50000; //N for integral
	int jmax = 10000; //iterations per radius
	int kmax = 40; //number of radii
	double boxsize = 30;
	double xval, yval, radius, radx, rady;
	double integral, outegral;;
	double norm=pow(sigma*sqrt(2*pi),1);

	ofstream r1;
	r1.open("1dradij.dat");


	for(int k=20;k<kmax;k++){
		radius = k*1.0;
		outegral=0;
		for(int j=0;j<jmax;j++){
			integral=0;
			radx = randDouble(-1*boxsize,boxsize);
			//rady = randDouble(-1*boxsize,boxsize);
			for(int i=0;i<imax;i++){
				xval = randDouble(-1*boxsize,boxsize);
				//yval = randDouble(-1*boxsize,boxsize);
				//integral+=gaussian(xval,sigma)*gaussian(yval,sigma);
				if(distcheck(xval,0,radx,0,radius)){
					integral +=gaussian(xval,sigma);
				}
			}
			outegral += gaussian(radx,sigma)*pow(2*boxsize,1)*integral/(norm*imax);
		}
		cout << pow(2*boxsize,1)*outegral/(norm*jmax) << endl;

		r1 << radius << " " << pow(2*boxsize,1)*outegral/(norm*jmax) << endl;
	}

	cout << pow(40,1)*integral/(norm*imax) << endl;

	r1.close();

}