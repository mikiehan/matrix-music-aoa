/*
 * Copyright 2016 <Admobilize>
 * MATRIX Labs  [http://creator.matrix.one]
 * This file is part of MATRIX Creator HAL
 *
 * MATRIX Creator HAL is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <map>
#include <string>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/FFT>

#include "./cross_correlation.h"
#include "./microphone_array.h"
#include "cpp/driver/microphone_array_location.h"

#define PI 3.14159265358979323846
//using namespace std;
//using namespace Eigen;

using namespace std;
using namespace Eigen;

typedef Matrix<float, 8, 128> Matrix8by128f;
typedef Matrix<float, 128, 8> Matrix128by8f;
typedef Matrix<float, 128, 128> Matrix128by128f;
typedef Matrix<float, 8, 8> Matrix8by8f;

typedef Matrix<complex<float>, 8, 8> MatrixC8by8f;
typedef Matrix<complex<float>, 1, 8> MatrixC1by8f;
typedef Matrix<MatrixC8by8f, 128, 1> MatrixCC128by1f;

typedef Matrix<float, 8, 7> Matrix8by7f;
typedef Matrix<float, 7, 8> Matrix7by8f;

typedef Matrix<complex<float>, 8, 7> MatrixC8by7f;
typedef Matrix<complex<float>, 7, 8> MatrixC7by8f;
typedef Matrix<float, 128, 90> Matrix128by90f;
typedef Matrix<float, 128, 38> Matrix128by38f;
typedef Matrix<float, 8, 1> Matrix8by1f;
typedef Matrix<complex<float>, 8, 1> MatrixC8by1f;
typedef Matrix<complex<float>, 8, 5> MatrixC8by5f;
typedef Matrix<float, 1, 1> Matrix1by1f;
typedef Matrix<complex<float>, 1, 1> MatrixC1by1f;
typedef Matrix<Matrix8by1f, 360, 1> MatrixN360by1f;
typedef Matrix<float, 360, 1> Matrix360by1f;
typedef Matrix<MatrixC8by1f, 360, 1> MatrixC360by1f;
typedef Matrix<MatrixC8by1f, 1, 360> MatrixC1by360f;

#include "cpp/driver/direction_of_arrival.h"



MatrixC360by1f steeringVector_;
//Matrix360by1f steeringVector_;

//typedef Matrix<float, 8, 128> Matrix8f;
//typedef Matrix<float, 8, 8> Matrixaf;

namespace matrix_hal {

DirectionOfArrival::DirectionOfArrival(MicrophoneArray& mics)
    : mics_(mics),
      length_(mics_.NumberOfSamples()),
      corr_(length_),
      current_mag_(4),
      current_index_(4),
      buffer_1D_(mics_.Channels() * mics_.SamplingRate()),
      buffer_2D_(mics_.Channels()),
      mic_direction_(0),
      azimutal_angle_(0.0),
      polar_angle_(0.0),
      snap_shots_(125.0), //1000 ms of data
      file_count_(0.0) {
  for (uint16_t c = 0; c < mics_.Channels(); c++) {
    buffer_2D_[c] = &buffer_1D_[c * mics_.SamplingRate()];
  }
}

void DirectionOfArrival::steeringVectorCalculate(){

    //using namespace std::complex_literals;
    //std::cout << std::fixed << std::setprecision(1);

    const std::complex<double> l(0, 1);

    for(int i = 0 ; i < 360 ; i ++){
	//cout << i << endl; 
    	MatrixC8by1f steeringVec_;
    	//steeringVec_(0,0) = 1.0;
    
	//steeringVector_(0,0) = 1;

    	for(int j = 0 ; j < 8 ; j++){
         
       		//steeringVec_(j,0) = std::exp(-l*(2*PI*16000*(0.054*(cos( (i*PI)/180.0) - cos( j*((360.0*PI)/(8.0*180.0)) - (i*PI)/180.0)))/334.0));
       		//steeringVec_(j,0) = std::exp( -l*( 2*PI*16000*( 0.054*(cos( (i*PI)/180.0 - (2*PI/8.0) ))/334.0)));
       		steeringVec_(j,0) = std::exp(l*( ((2*PI*16000*0.054)/335.0)*cos( (i*PI/180.0- (2*PI*PI*j)/(8.0*180.0)) ) ) );
       		//steeringVec_(j,0) = std::exp(-l*(2*PI*16000*(0.054*(cos( (i*PI)/180.0) - cos( j*((360.0*PI)/(8.0*180.0))))/334.0)));
		//steeringVector_(i,j) = cos(2*PI*16000*((0.054*(cos(i*PI/180.0) - cos( j*(360.0*PI/8.0*180.0) - i*PI/180.0)))/330.0));
    	}

	steeringVector_(i,0) = MatrixC8by1f(steeringVec_);
     }  
}

void DirectionOfArrival::Calculate() {


  //int max_tof = mics_.NumberOfSamples();
  //int max_tof = 10; //Ignore thes initial samples

 //cout << mics_.NumberOfSamples() << endl;
 //cout << mics_.Channels() << endl;
  //cout << mics_.NumberOfSamples() << endl;

  //Create the Matrix
  MatrixC8by8f musicFinalMatrix; 
  MatrixCC128by1f musicFFTEachBinMatrix; 

  for(int i = 0; i < snap_shots_; i++){

        FFT<float> fft;
  	Matrix<float, 8, 128> musicMatrix;
	Matrix<complex<float>, 8, 128> outMatrix;
        outMatrix.setZero(8, 128);
      

  	//Putting the raw Data into buffers
  	for (uint32_t s = 0; s < mics_.NumberOfSamples(); s++) {
    	
		for (uint16_t c = 0; c < mics_.Channels(); c++) { /* mics_.Channels()=8 */
      			//buffer_2D_[c][s] = mics_.At(s, c); //Delayed data
      			buffer_2D_[c][s] = mics_.Atr(s, c); //Raw data
      			//musicMatrix(s,c) = mics_.Atr(s, c); //Raw data
    		}

  	}
        
  	for (uint32_t c = 0; c < mics_.Channels(); c++) {
    		for (uint16_t s = 0; s < mics_.NumberOfSamples(); s++) { // mics_.Channels()=8
			musicMatrix(c,s) = buffer_2D_[c][s];
		}	
  	}

        //FFT Operations


	for (int k = 0; k < musicMatrix.rows(); k++) {
    		Matrix<complex<float>, 1, 128> tmpOut;
    		fft.fwd(tmpOut, musicMatrix.row(k));
    		outMatrix.row(k) = tmpOut;
	}

	for (int k = 0; k < musicMatrix.cols(); k++) {
    		Matrix<complex<float>, 8, 1> tmpOut;
    		fft.fwd(tmpOut, outMatrix.col(k));
    		outMatrix.col(k) = tmpOut;
	}

        //Reading again Mic Values
  	mics_.Read();

        //Calculating the avearge matrix values
        if(i == 0 ){

                Matrix<complex<float>, 128, 8> outMatrixAdj(outMatrix.adjoint());
  		musicFinalMatrix = outMatrix*(outMatrixAdj);
  		//musicFinalMatrix = outMatrix*(outMatrix.adjoint());
        
               /*
		for(int z=0; z < outMatrix.cols(); z++){
		
			MatrixC8by1f tmpMatrix;
			tmpMatrix = outMatrix.col(z);
                        MatrixC8by8f tmpMatrix1;
                        MatrixC1by8f tmpMatrix2(tmpMatrix.adjoint());
                        tmpMatrix1 = tmpMatrix*tmpMatrix2;
			musicFFTEachBinMatrix(z,0) = tmpMatrix1;
                }*/

	}else{

		MatrixC8by8f tmpMtx;
                Matrix<complex<float>, 128, 8> outMatrixAdj(outMatrix.adjoint());
  		tmpMtx = musicFinalMatrix + (outMatrix*(outMatrixAdj));
		musicFinalMatrix = tmpMtx;
		
		//musicFinalMatrix = musicFinalMatrix + (outMatrix*(outMatrix.adjoint()));

                /*
		for(int z=0; z < outMatrix.cols(); z++){
		
			MatrixC8by1f tmpMatrix;
			tmpMatrix = outMatrix.col(z);
			MatrixC1by8f tmpMatrix1(tmpMatrix.adjoint());
                        MatrixC8by8f othMatrix1(tmpMatrix*tmpMatrix1);
                        
			MatrixC8by8f othMatrix;
			othMatrix = ( musicFFTEachBinMatrix(z,0) + othMatrix1 );
                        musicFFTEachBinMatrix(z,0) = othMatrix;
                }*/
	}
   
	//cout << "----------------------" << endl << musicFinalMatrix << endl << "----------------------" << endl;
   }

   //musicFinalMatrix.setZer(8,8);
   musicFinalMatrix = (1.0/(1.0*snap_shots_))*musicFinalMatrix;

 /*  for(int y = 0; y < musicFFTEachBinMatrix.rows(); y++ ){

	MatrixC8by8f tmpMatrix;
        tmpMatrix = musicFFTEachBinMatrix(y,0);
	MatrixC8by8f tmpMatrix1;
	tmpMatrix1 =  (1.0/(1.0*snap_shots_))*tmpMatrix;
	musicFFTEachBinMatrix(y,0) = tmpMatrix1;	
	
   }*/

   //cout << "F----------------------" << endl << musicFinalMatrix << endl << "----------------------" << endl;
   
   //For general Mtarix
   //EigenSolver<MatrixC8by8f> es(musicFinalMatrix);
   
   //For real symmetric matrix
   //SelfAdjointEigenSolver<MatrixC8by8f> es(musicFinalMatrix);
  
   //For Complex Mtarix
   ComplexEigenSolver<MatrixC8by8f> es;
   es.compute(musicFinalMatrix);

   
   cout << "The eigenvalues of Correlation Matrix are:" << endl << es.eigenvalues() << endl;
   
   /*
    ofstream myfile;
  
   myfile.open ("example_" + std::to_string(file_count_) + ".dat");
   for(int y = 0; y < musicFFTEachBinMatrix.rows(); y++ ){

	MatrixC8by8f tmpMatrix;
	tmpMatrix = musicFFTEachBinMatrix(y,0);

        ComplexEigenSolver<MatrixC8by8f> esTemp;
        esTemp.compute(tmpMatrix);

        MatrixC8by1f egVal = esTemp.eigenvalues();
        
        std::string eigLine("");

        for(int a = 0; a < 8; a++){

		double val = abs(egVal(a,0));
                
                std::string vstr = std::to_string(log10(val));
                if( a == 0 ){
                   eigLine = vstr;
                }else{
                      eigLine = eigLine + " " + vstr;
                }
        }
   	myfile << y << " " << eigLine << endl;
	//cout << endl << "------The eigenvalues of Correlation Matrix " << y << "th column are : " << endl << esTemp.eigenvalues() << endl;
        //cout << endl << "---------------------------------------------------------------" << endl;
   }
        myfile.close();

   */

   MatrixC8by7f Qmat;

   for(int l=0; l<7 ; l++)
      Qmat.col(l) = es.eigenvectors().col(l);
 
   //cout << "The matrix of eigenvectors, Q, is:" << endl << Qmat << endl << endl;

  //VectorXd m_solved_val = es.eigenvalues.real(); 

   Matrix360by1f phaseSpectrum_;

  //ofstream myfile;
  
  //myfile.open ("example_" + std::to_string(file_count_) + ".dat");
 
  MatrixC7by8f QmatAdj(Qmat.adjoint());

  for(int i = 0; i < 360; i++){

   MatrixC1by8f stVecAdj(steeringVector_(i,0).adjoint());

   MatrixC1by1f val = stVecAdj*Qmat*QmatAdj*steeringVector_(i,0);
   //MatrixC1by1f val = steeringVector_(i,0).adjoint()*Qmat*Qmat.adjoint()*steeringVector_(i,0);

   phaseSpectrum_(i,0) = 1.0/(val.norm());
   //myfile << i << " " << 1.0/(val.norm()) << endl;
 }

  //myfile.close();
  file_count_ += 1;
  
  //get location of maximum
  MatrixXf::Index maxRow, maxCol;
  float max = phaseSpectrum_.maxCoeff(&maxRow, &maxCol);

  //cout << phaseSpectrum_ << endl;

    cout << "Value :" << max << "; Max Degree : " << maxRow << endl;

   if ( max > 0.05) {
  
    mic_direction_ = maxRow/45;
    azimutal_angle_ = atan2(micarray_location[mic_direction_][1],
                            micarray_location[mic_direction_][0]);
    //polar_angle_ = fabs(index) * M_PI / 2.0 / float(max_tof - 1);

  }

  //Sorting Eigen values
  //VectorXf eVecs = es.eigenvectors().real(); // Keep only the real part of complex matrix
  //VectorXf eVals = es.eigenvalues().real(); // Keep only the real part of complex matrix

  // Sort by descending eigenvalues:
   /*
  std::vector<std::pair<Scalar,Index> > D;
  D.reserve(eVals.size());

  for (Index i=0;i<eVals.size();i++)
	D.push_back(std::make_pair<Scalar,Index>(eVals.coeff(i,0),i));

   std::sort(D.begin(),D.end());

   MATRIX1 sortedEigs;

   sortedEigs.resizeLike(eVecs);

   for (int i=0;i<eVals.size();i++)
   {
          eVals.coeffRef(i,0)=D[i].first;
          sortedEigs.col(i)=eVecs.col(D[i].second);
   }

     eVecs = sortedEigs;
*/

   //FullPivLU<Matrix8by8f> lu_decomp(musicFinalMatrix);
   //auto rank = lu_decomp.rank();
   //cout << "Rank " << rank << endl; 
    
  /*

  //Finding the correlation between every 4th i.e. opposite mics

  for (int channel = 0; channel < 4; channel++) {
    corr_.Exec(buffer_2D_[channel + 4], buffer_2D_[channel]);

    float* c = corr_.Result();

    //for (int i = 0; i< mics_.NumberOfSamples(); i++ )   
    //	cout << c[i] <<endl;
    //cout << endl;
  //}

    int index = max_tof;
    float m = c[index];

    //Find the max correlation and its corresponding index
    for (int i = 1; i < length_; i++)
      if (c[i] > m) {
        index = i;
        m = c[i];
      }

    for (int i = max_tof; i < length_; i++)
      if (c[i] > m) {
        index = i;
        m = c[i];
      }

    current_mag_[channel] = m;
    current_index_[channel] = index;
  }

  for( int x = 0; x < 4 ; x++){
    cout << " Magnitude : " << current_mag_[x] << " Index : " << current_index_[x] << endl; 
  }
 
  int dir = 0;
  int index = current_index_[0];
  float mag = current_mag_[0];
  for (int channel = 1; channel < 4; channel++) {
    if (mag < current_mag_[channel]) {
      dir = channel;
      mag = current_mag_[channel];
      index = current_index_[channel];
    }
  }

 if (index > 64) index = -(128 - index);
  
  std::cout << "Min Mag = " << mag << "Min Index "<< index <<"\n";

if (mag > 3.8e8) {
    if (index < 0) {
      mic_direction_ = (dir + 4);
    } else {
      mic_direction_ = dir;
    }

    azimutal_angle_ = atan2(micarray_location[mic_direction_][1],
                            micarray_location[mic_direction_][0]);
    polar_angle_ = fabs(index) * M_PI / 2.0 / float(max_tof - 1);

  }

 */

}

};  // namespace matrix_hal
