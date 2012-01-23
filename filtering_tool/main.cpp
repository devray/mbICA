#include <armadillo>
#include <functional>
#include <fftw3.h>
#include <iostream>
#include <string>

using namespace arma;
using namespace std;

class Transformer {
public:
	Transformer(int N);
	~Transformer();
	void filterSignal(const mat &sig, int column, 
					  int startFreq, int stopFreq, mat &out);
private:
	int getSampleNumber(int freq);
	int N_;
	fftw_complex *transBuffer_;
	fftw_complex *filterBuffer_;
	double *dataBuffer_;
	double *innerBuffer_;
	fftw_plan forwardPlan_;
	fftw_plan backwardPlan_;
};

int main(int argc, char *argv[]) {
	if(argc < 4) {
		cerr << "Too few parameters" << endl;
		return -1;
	}
	mat signal;
	string fileName = argv[1];
	cout << "Opening " << fileName << "..." << endl;
    signal.load(argv[1], raw_ascii);
	cout << "Opened "  << signal.n_cols << " signals with " << signal.n_rows << " samples each" << endl;

	int startFreq = atoi(argv[2]);
	int stopFreq = atoi(argv[3]);

	cout << "Filtering from " << startFreq << "Hz to " << stopFreq << "Hz." << endl;
	Transformer t(signal.n_rows);

	for(unsigned i = 0; i < signal.n_cols; ++i) {
		cout << "\tFiltering " << i << "th signal." << endl;
		t.filterSignal(signal, i, startFreq, stopFreq, signal);
	}

	cout << "Saving..." << endl;
	signal.save(("filtered_" + fileName).c_str(), raw_ascii);
	return 0;
}

Transformer::Transformer(int N)
	: N_(N) {
	// Creating buffers
	transBuffer_ = (fftw_complex *) fftw_malloc(N_ * sizeof(fftw_complex));
	filterBuffer_ = (fftw_complex *) fftw_malloc(N_ * sizeof(fftw_complex));

	dataBuffer_ = (double *) fftw_malloc(N_ * sizeof(double));
	innerBuffer_ = (double *) fftw_malloc(N_ * sizeof(double));

	// Creating plans
	forwardPlan_ = fftw_plan_dft_r2c_1d(N_, dataBuffer_, transBuffer_, FFTW_MEASURE);
	backwardPlan_ = fftw_plan_dft_c2r_1d(N_, filterBuffer_, innerBuffer_, FFTW_MEASURE);
}

Transformer::~Transformer() {
	// Destorying plans
	fftw_destroy_plan(forwardPlan_);
	fftw_destroy_plan(backwardPlan_);

	// Destorying buffers
	fftw_free(transBuffer_);
	fftw_free(filterBuffer_);
	fftw_free(dataBuffer_);
	fftw_free(innerBuffer_);
}

void Transformer::filterSignal(const mat &sig, int column, 
							   int startFreq, int stopFreq, mat &out) {
	int i;
	// copying from input
	for(i = 0; i < N_; ++i) {
		dataBuffer_[i] = sig(column, i);
	}

	// converting freq to sample nnumbers
	startFreq = getSampleNumber(startFreq);
	stopFreq = getSampleNumber(stopFreq);

	// FFT
	fftw_execute(forwardPlan_);

	// Filtering - bandpass
	for(i = 0; i < startFreq; ++i) {
		filterBuffer_[i][0] = filterBuffer_[i][1] = 0.0;
	}
	for(; i < stopFreq; ++i) {
		filterBuffer_[i][0] = transBuffer_[i][0];
		filterBuffer_[i][1] = transBuffer_[i][1];
	}
	for(; i < N_; ++i) {
		filterBuffer_[i][0] = filterBuffer_[i][1] = 0.0;
	}

	// inverse FFT
	fftw_execute(backwardPlan_);

	// copying to output
	for(i = 0; i < N_; ++i) {
		out(column, i) = innerBuffer_[i] / N_;
	}
}

int Transformer::getSampleNumber(int freq) {
	return ((freq * N_) >> 7);
}

