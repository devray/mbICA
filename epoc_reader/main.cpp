#include <libepoc.h>
#include <iostream>
#include <armadillo>
#include <signal.h>

using namespace std;
using namespace arma;

bool stop = false;

void set_stop(int) {
	stop = true;
}

int main(int argc, char *argv[]) {
	unsigned int row_num = 0;

	if(argc > 1) {
		row_num = atoi(argv[1]);
		cout << "Recording " << argv[1] << " samples" << endl;
	}

	signal(SIGINT, &set_stop);

	if(epoc_get_count(EPOC_VID, EPOC_PID) == 0){
		cerr << "Cannot find device with vid: " << EPOC_VID << "; pid: " << EPOC_PID << endl;
		return 0;
	}

	epoc_device *device = epoc_open(EPOC_VID, EPOC_PID, 0);
	if(device == NULL) {
		cerr << "Cannot open the device" << endl;
		return 0;
	}
	cout << "EEGDevice opened..." << endl;

	epoc_handler *handler;
	if((handler = epoc_init(device, RESEARCH_HEADSET)) == NULL) {
		cerr << "Cannot init the device!" << endl;
		epoc_close(device);
		return 0;
	}

	mat data(row_num, 16);
	
	struct epoc_frame fr;
	for(unsigned ix = 0; !stop && (row_num == 0 || ix < row_num); ++ix) {
		epoc_get_next_frame(handler, &fr);

		if(row_num == 0) {
			data.resize(ix + 1, 16);
		}

		for(int i=0; i < 16; ++i) {
			data(ix, i) = fr.electrode[i];
		}
	}

	data.save("data.dat", arma_ascii);
	data.save("data.octave", raw_ascii);

	return 0;
}
