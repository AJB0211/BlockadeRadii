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





int main(int argc, char * const argv[]){
	int imax = 50000; //N for integral
	int jmax = 10000; //iterations per radius
	int kmax = 1; //number of radii
	double boxsize = 30;
	double xval, yval, radius, radx, rady;
	double integral, outegral;;
	double norm=pow(sigma*sqrt(2*pi),1);

	//ofstream r1;
	//r1.open("1dradij.dat");


	for(int k=0;k<kmax;k++){
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

		//r1 << radius << " " << pow(2*boxsize,1)*outegral/(norm*jmax) << endl;
	}

	cout << pow(40,1)*integral/(norm*imax) << endl;

	//r1.close();

}