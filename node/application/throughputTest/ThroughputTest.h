/****************************************************************************
 *  Copyright: National ICT Australia,  2007 - 2010                         *
 *  Developed at the ATP lab, Networked Systems research theme              *
 *  Author(s): Athanassios Boulis, Yuriy Tselishchev                        *
 *  This file is distributed under the terms in the attached LICENSE file.  *
 *  If you do not find this file, copies can be found by writing to:        *
 *                                                                          *
 *      NICTA, Locked Bag 9013, Alexandria, NSW 1435, Australia             *
 *      Attention:  License Inquiry.                                        *
 *                                                                          *
 ****************************************************************************/
#ifndef _THROUGHPUTTEST_H_
#define _THROUGHPUTTEST_H_

#include "VirtualApplication.h"
#include <map>
#include <iostream>
#include<cmath>
#include<algorithm>

#include<ctime>     //时间的创建
#include<cassert>   //定义
using namespace std;

enum ThroughputTestTimers {
	SEND_PACKET = 1,
	ADD_PACKET = 2,
	ENERGY_TRAN = 3,
	ENERGY_COL = 4,
	SHOW_ENERGY = 5,
	SINK_SEND = 6,
	WRITE_DATA = 7,
	CHANGE_SLOTS = 8,
	TEST_SLOTS = 9
};

class ThroughputTest: public VirtualApplication {
 public:	
	double packet_rate;
	double startupDelay;
	double delayLimit;
	float packet_spacing;
	int dataSN;
	int recipientId;
	string recipientAddress;
	int packet_cap = 50000;
	int packet_count=0;
	float send_spacing;
	int send_rate=4;
	int max_energy_tx;
	//variables below are used to determine the packet delivery rates.	
	int numNodes;
	double node_two_consume=0;
	map<long,int> packetsReceived;
	map<long,int> bytesReceived;
	map<long,int> packetsSent;
	vector<vector<double>> entropyData;
	vector<vector<double>> entropyData2;

	vector<double> w={0.33333,0.33333,0.33333};
	vector<double> w2={0.33333,0.33333,0.33333};

	vector<double> energy_harvest={49008.747,46692.35,63286.47,51347.867,44753.79,29372.405,55129.27,17308.6,21942.99,41480.61};
	int energy_index=0;
	int energy_index2=0;
	int packageCapacity=245;
	double money=0;
	double addMoney=100;
	int sendCount=0;
	//nJ nJ/bit!!
	double E_TX=36.1;
	double E_RX=16.7;
	double E_FS=0.01;
	//10*10^-12 in J,10*10^-3 in nJ;
	double E_MP=0.0000013;
	double dmp=87.706;
	vector<int> slot_count={3,3,2,2,2,2,1,1};
	int base_send_rate=50;
	
 protected:
	void startup();
	void fromNetworkLayer(ApplicationPacket *, const char *, double, double);
	void handleRadioControlMessage(RadioControlMessage *);
	void timerFiredCallback(int);
	void finishSpecific();
	

 public:
	int getPacketsSent(int addr) { return packetsSent[addr]; }
	int getPacketsReceived(int addr) { return packetsReceived[addr]; }
	int getBytesReceived(int addr) { return bytesReceived[addr]; }
	
};

#endif				// _THROUGHPUTTEST_APPLICATIONMODULE_H_
