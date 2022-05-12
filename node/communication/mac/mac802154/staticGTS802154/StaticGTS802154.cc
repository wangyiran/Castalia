/****************************************************************************
 *  Copyright: National ICT Australia,  2007 - 2012                         *
 *  Developed at the ATP lab, Networked Systems research theme              *
 *  Author(s): Yuriy Tselishchev                                            *
 *  This file is distributed under the terms in the attached LICENSE file.  *
 *  If you do not find this file, copies can be found by writing to:        *
 *                                                                          *
 *      NICTA, Locked Bag 9013, Alexandria, NSW 1435, Australia             *
 *      Attention:  License Inquiry.                                        *
 *                                                                          *
 ****************************************************************************/

#include "StaticGTS802154.h"

Define_Module(StaticGTS802154);

/***
 * Initialising some parameters, specific to Static GTS module
 * by overriding the startup() method. Important to call startup()
 * of the parent module in the end, otherwise it will not initialize
 ***/
void StaticGTS802154::startup() {
	// initialise GTS-specific parameters
	GTSlist.clear(); totalGTS = 0; assignedGTS = 0;
	requestGTS = par("requestGTS");
	gtsOnly = par("gtsOnly");

	// other parameters are from Basic802154, need to read them for GTS scheduling
	totalSlots = par("numSuperframeSlots"); 	
	baseSlot = par("baseSlotDuration");
	minCap = par("minCAPLength");
	frameOrder = par("frameOrder");
	return Basic802154::startup();
	setTimer(SHOW_GTS,10);
}


void StaticGTS802154::timerFiredCallback(int index)
{
	switch (index) {
		
		// Start of a new superframe
		case FRAME_START: {
			
			GTSlist.clear(); totalGTS = 0;
			GTSlist.resize(7);
			ofstream outfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/wangzha.txt", ios::app);

			outfile<<" frame show gts list : "<<totalGTS<<endl;
	

			ifstream myfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/gts"); 
			string slotsString;
			getline(myfile,slotsString);
			outfile<<slotsString<<endl;
			vector<int> slt(8,0);
			int idx=0;
			for(int k=0;k<slotsString.size();k++){
				if(slotsString[k]!=' '){
					slt[idx]=(slotsString[k]-'0');
					idx++;
				}
			}

			int start=slt[0];
			for (int i = 0; i < (int)GTSlist.size(); i++) {
				GTSlist[i].owner=i+1;
				GTSlist[i].length=slt[i+1];
				GTSlist[i].start=start;
				totalGTS+=slt[i+1];
				start+=slt[i+1];
				if(self==i+1){
					assignedGTS=slt[i+1];
				}
			}

			vector<Basic802154GTSspec>::iterator i;
			for (i = GTSlist.begin(); i != GTSlist.end(); i++) {
				outfile<<"slot owner:"<<i->owner<<" start:"<<i->start<<" len:"<<i->length<<"meishuchu"<<endl;
			}

			//myfile.close();
			outfile.close();

			if (isPANCoordinator) {	// as a PAN coordinator, create and broadcast beacon packet
				beaconPacket = new Basic802154Packet("PAN beacon packet", MAC_LAYER_PACKET);
				beaconPacket->setDstID(BROADCAST_MAC_ADDRESS);
				beaconPacket->setPANid(SELF_MAC_ADDRESS);
				beaconPacket->setMac802154PacketType(MAC_802154_BEACON_PACKET);
				beaconPacket->setBeaconOrder(beaconOrder);
				beaconPacket->setFrameOrder(frameOrder);
				if (++macBSN > 255) macBSN = 0;
				beaconPacket->setBSN(macBSN);
				beaconPacket->setCAPlength(numSuperframeSlots);
				
				// GTS fields and CAP length are set in the decision layer
				prepareBeacon_hub(beaconPacket);
				
				beaconPacket->setByteLength(BASE_BEACON_PKT_SIZE + beaconPacket->getGTSlistArraySize() * GTS_SPEC_FIELD_SIZE);
				CAPlength = beaconPacket->getCAPlength();
				CAPend = CAPlength * baseSlotDuration * (1 << frameOrder) * symbolLen;
				sentBeacons++;

				trace() << "Transmitting [PAN beacon packet] now, BSN = " << macBSN;
				setMacState(MAC_STATE_CAP);
				toRadioLayer(beaconPacket);
				toRadioLayer(createRadioCommand(SET_STATE, TX));
				setTimer(ATTEMPT_TX, TX_TIME(beaconPacket->getByteLength()));
				beaconPacket = NULL;

				currentFrameStart = getClock() + phyDelayRx2Tx;
				setTimer(FRAME_START, beaconInterval * symbolLen);
			} else {	// if not a PAN coordinator, then wait for beacon
				toRadioLayer(createRadioCommand(SET_STATE, RX));
				setTimer(BEACON_TIMEOUT, guardTime * 3);
			}
			break;
		}

		case GTS_START: {

		
			if (macState == MAC_STATE_SLEEP) {
				toRadioLayer(createRadioCommand(SET_STATE, RX));
			}
			setMacState(MAC_STATE_GTS);
			
			// we delay transmission attempt by the time requred by radio to wake up
			// note that GTS_START timer was scheduled exactly phyDelaySleep2Tx seconds
			// earlier than the actual start time of GTS slot
			setTimer(ATTEMPT_TX, phyDelaySleep2Tx);

			// set a timer to go to sleep after this GTS slot ends
			setTimer(SLEEP_START, phyDelaySleep2Tx + GTSlength);

			// inform the decision layer that GTS has started
			startedGTS_node();
			break;
		}

		// beacon timeout fired - indicates that beacon was missed by this node
		case BEACON_TIMEOUT: {
			lostBeacons++;
			if (lostBeacons >= maxLostBeacons) {
				trace() << "Lost synchronisation with PAN " << associatedPAN;
				setMacState(MAC_STATE_SETUP);
				associatedPAN = -1;
				desyncTimeStart = getClock();
				disconnectedFromPAN_node();
				if (currentPacket) clearCurrentPacket("No PAN");
			} else if (associatedPAN != -1) {
				trace() << "Missed beacon from PAN " << associatedPAN <<
				    ", will wake up to receive next beacon in " <<
				    beaconInterval * symbolLen - guardTime * 3 << " seconds";
				setMacState(MAC_STATE_SLEEP);
				toRadioLayer(createRadioCommand(SET_STATE, SLEEP));
				setTimer(FRAME_START, beaconInterval * symbolLen - guardTime * 3);
			}
			break;
		}
		
		// packet was not received
		case ACK_TIMEOUT: {
			collectPacketHistory("NoAck");
			attemptTransmission("ACK timeout");
			break;
		}

		// previous transmission is reset, attempt a new transmission
		case ATTEMPT_TX: {
			if (getTimer(ACK_TIMEOUT) != -1) break;
			attemptTransmission("ATTEMPT_TX timer");
			break;
		}

		// timer to preform Clear Channel Assessment (CCA) 
		case PERFORM_CCA: {
			if (macState == MAC_STATE_GTS || macState == MAC_STATE_SLEEP) break;
			CCA_result CCAcode = radioModule->isChannelClear();
			if (CCAcode == CLEAR) {
				//Channel clear
				if (--CW > 0) {
					setTimer(PERFORM_CCA, unitBackoffPeriod * symbolLen);
				} else {
					transmitCurrentPacket();
				}
			} else if (CCAcode == BUSY) {
				//Channel busy
				CW = enableSlottedCSMA ? 2 : 1;
				if (++BE > macMaxBE)
					BE = macMaxBE;
				if (++NB > macMaxCSMABackoffs) {
					collectPacketHistory("CSfail");
					currentPacketRetries--;
					attemptTransmission("Current NB exeeded maxCSMAbackoffs");
				} else {
					performCSMACA();
				}
			} else if (CCAcode == CS_NOT_VALID_YET) {
				//Clear Channel Assesment (CCA) pin is not valid yet
				setTimer(PERFORM_CCA, phyDelayForValidCS);
			} else {	
				//Clear Channel Assesment (CCA) pin is not valid at all (radio is sleeping?)
				trace() << "ERROR: isChannelClear() called when radio is not ready";
				toRadioLayer(createRadioCommand(SET_STATE, RX));
			}
			break;
		}

		case SLEEP_START: {
			// SLEEP_START timer can sometimes be scheduled in the end of a frame
			// i.e. when BEACON_ORDER = FRAME_ORDER, overlapping with the interval 
			// when a node already tries to prepare for beacon reception. Thus 
			// check if BEACON_TIMEOUT timer is set before going to sleep
			if (getTimer(BEACON_TIMEOUT) != -1) break;

			cancelTimer(PERFORM_CCA);
			setMacState(MAC_STATE_SLEEP);
			toRadioLayer(createRadioCommand(SET_STATE, SLEEP));
			break;
		}
		
		case BACK_TO_SETUP: {
			// This timer is scheduled to the end of the CAP period 
			// when beacon is received, but node is not (yet) connected.
			// So when this timer fires and node is not connected, it 
			// has to go back to setup stage
			if (associatedPAN == -1) setMacState(MAC_STATE_SETUP);
		}
		
		case SHOW_GTS:{
			for (int i = 0; i < (int)GTSlist.size(); i++) {
				GTSlist[i].length=1;

			}

		}
		
	}
}


/***
 * GTS request received by hub, need to return the number of 
 * slots to be granted. Can return 0 to indicate request denial
 ***/
int StaticGTS802154::gtsRequest_hub(Basic802154Packet *gtsPkt) {
	// //Length of CAP after lengths of all GTS slots are subtracted
	// int CAPlength = totalSlots;
	// 			ofstream outfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/wangzha.txt", ios::app);

	// 		outfile<<"func"<<GTSlist.size()<<endl;
	// outfile.close();

	// //check if the node already exists in the GTS list
	// vector<Basic802154GTSspec>::iterator i;
	// int total = 0;
	// for (i = GTSlist.begin(); i != GTSlist.end(); i++) {
	// 	total++;
	// 	if (i->owner == gtsPkt->getSrcID()) {
	// 		if (i->length == requestGTS) {
	// 			return i->length;
	// 		} else {
	// 			totalGTS -= i->length;
	// 			GTSlist.erase(i);
	// 			total--;
	// 		}
	// 	} else {
	// 		CAPlength -= i->length;
	// 	}
	// }
	
	// //node not found, or requested slots changed
	// if (total >= 7 || (CAPlength - gtsPkt->getGTSlength()) *
	//     baseSlot * (1 << frameOrder) < minCap) {
	// 	trace() << "GTS request from " << gtsPkt->getSrcID() <<
	// 	    " cannot be acocmodated";
	// 	return 0;
	// }
	
	// Basic802154GTSspec newGTSspec;
	// newGTSspec.length = gtsPkt->getGTSlength();
	// totalGTS += newGTSspec.length;
	// newGTSspec.owner = gtsPkt->getSrcID();
	// GTSlist.push_back(newGTSspec);
	// outfile<<"len:"<<totalGTS<<endl;
	// TT=GTSlist;
	// outfile<<"TTlen:"<<TT.size()<<endl;
	// 	outfile.close();

	// return newGTSspec.length;
}

/***
 * Hub can alter the beacon before broadcasting it
 * In particular, assign GTS slots and set CAP length
 ***/
void StaticGTS802154::prepareBeacon_hub(Basic802154Packet *beaconPacket) {
	int CAPlength = totalSlots;
	beaconPacket->setGTSlistArraySize(GTSlist.size());
	for (int i = 0; i < (int)GTSlist.size(); i++) {
		if (CAPlength > GTSlist[i].length) {
			CAPlength -= GTSlist[i].length;
			GTSlist[i].start = CAPlength + 1;
			beaconPacket->setGTSlist(i, GTSlist[i]);
		} else {
			trace() << "Internal ERROR: GTS list corrupted";
			GTSlist.clear(); totalGTS = 0;
			beaconPacket->setGTSlistArraySize(0);	
			CAPlength = totalSlots;
			break;
		}
	}
	beaconPacket->setCAPlength(CAPlength);
}

/***
 * If disconnected from PAN, also need to reset GTS slots
 ***/
void StaticGTS802154::disconnectedFromPAN_node() {
	assignedGTS = 0;
}

/***
 * GTS request was successful
 ***/
void StaticGTS802154::assignedGTS_node(int slots) {
	assignedGTS = slots;
}

/***
 * Transmission of data packet requested earlier is complete
 * status string holds comma separated list of outcomes
 * for each transmission attempt
 ***/
void StaticGTS802154::transmissionOutcome(Basic802154Packet *pkt, bool success, string status) {
	if (getAssociatedPAN() != -1) {
		if (assignedGTS == 0 && requestGTS > 0) {
			transmitPacket(newGtsRequest(getAssociatedPAN(), requestGTS));
		} else if (TXBuffer.size()) {
			Basic802154Packet *packet = check_and_cast<Basic802154Packet*>(TXBuffer.front());
			TXBuffer.pop();
			transmitPacket(packet,0,gtsOnly);
		}
	}
}

bool StaticGTS802154::acceptNewPacket(Basic802154Packet *newPacket) 
{
	if (getAssociatedPAN() != -1 && getCurrentPacket() == NULL) {
		transmitPacket(newPacket,0,gtsOnly);
		return true;
	}
	return bufferPacket(newPacket);
}

/***
 * Timers can be accessed by overwriting timerFiredCallback
 **/
/*
void StaticGTS802154::timerFiredCallback(int index) {
	switch(index) {
		case NEW_TIMER: {
			//do something
			break;
		}
		
		default: {
			//important to call the function of the parent module
			Basic802154::timerFiredCallback(index);
		}
	}
}
*/

