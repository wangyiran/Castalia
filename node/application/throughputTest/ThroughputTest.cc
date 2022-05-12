/****************************************************************************
 *  Copyright: National ICT Australia,  2007 - 2011                         *
 *  Developed at the ATP lab, Networked Systems research theme              *
 *  Author(s): Athanassios Boulis, Yuriy Tselishchev                        *
 *  This file is distributed under the terms in the attached LICENSE file.  *
 *  If you do not find this file, copies can be found by writing to:        *
 *                                                                          *
 *      NICTA, Locked Bag 9013, Alexandria, NSW 1435, Australia             *
 *      Attention:  License Inquiry.                                        *
 *                                                                          *
 ****************************************************************************/

#include "ThroughputTest.h"

Define_Module(ThroughputTest);




vector<double> dimensionlessTreatment(vector<vector<double>> &inputData,vector<double> &eachMax,vector<double> &eachMin,vector<double> &eachDiff){
    vector<double> sum_ij(eachMax.size(),0);
    // + - +
    //Positive indicator
    for(int i=0;i<inputData.size();i++){
        double X_ij = (inputData[i][0] - eachMin[0])/eachDiff[0];
       // double X_ij = (eachMax[0] - inputData[i][0])/eachDiff[0];
        sum_ij[0]+=X_ij;
        inputData[i][0]=X_ij+0.00001;
    }

    //Negative indicator
    
	for(int i=0;i<inputData.size();i++){
		double X_ij = (eachMax[1]-inputData[i][1])/eachDiff[1];
		sum_ij[1]+=X_ij;
		inputData[i][1]=X_ij+0.00001;
	}

    //Positive indicator
	for(int i=0;i<inputData.size();i++){
		double X_ij = (eachMax[2]-inputData[i][2])/eachDiff[2];
		sum_ij[2]+=X_ij;
		inputData[i][2]=X_ij+0.00001;
	}
    return sum_ij;
}

void calculateProbability(vector<vector<double>> &inputData,vector<double> sum_ij){
    for(int i=0;i<3;i++){
        double this_sum=sum_ij[i]+0.0001;
        for(int j=0;j<inputData.size();j++){
            inputData[j][i]/=this_sum;
        }
    }
}

vector<double> calculateEntropy(vector<vector<double>> &inputData){
    vector<double> ej(inputData[0].size(),0);
    for(int i=0;i<3;i++){
        for(int j=0;j<inputData.size();j++){
            double pij=inputData[j][i];
			double temp=pij*log(pij);
			if (isnan(temp)){
				temp=0;
			}
            ej[i]+=temp;
        }
    }
    int n=inputData.size();
    for(int i=0;i<3;i++){
        ej[i]=-(ej[i]/log(n));
    }
    return ej;
}

vector<double> entropyMethod(vector<vector<double>> inputData){
    //Data preprocessing
        vector<double> eachMax,eachMin,eachDiff;
    int paramsCount=3;
    for(int i=0;i<paramsCount;i++){
        double tempMax=INT32_MIN;
        double tempMin=INT32_MAX;
        for(auto data :inputData){
            tempMax=data[i]>tempMax?data[i]:tempMax;
            tempMin=data[i]<tempMin?data[i]:tempMin;
        }
        eachMax.push_back(tempMax);
        eachMin.push_back(tempMin);
        eachDiff.push_back(tempMax-tempMin);
    
    }   

    //Dimensionless treatment
    vector<double> sum_ij(paramsCount,0);
    sum_ij=dimensionlessTreatment(inputData,eachMax,eachMin,eachDiff);

    //Calculate probability
    calculateProbability(inputData,sum_ij);

    //Calculate entropy
    vector<double> e_j(paramsCount,0);
    e_j=calculateEntropy(inputData);
	cout<<"e_j:"<<e_j[0]<<' '<<e_j[1]<<' '<<e_j[2]<<endl;

    //Calculate gj
    double sum_gj=0;
    for(int i=0;i<paramsCount;i++){
        e_j[i]=1-e_j[i];
        sum_gj+=e_j[i];
    }

    //Calculate weight
    vector<double> wj;
    for(int i=0;i<paramsCount;i++){
        wj.push_back(e_j[i]/sum_gj);
    }

    for(int i=0;i<paramsCount;i++){
        cout<<"weight"<<i<<"="<<wj[i]<<endl;
    }

    return wj;
}




double cal_dir(int s,int d)
{
	cTopology *topo;	// temp variable to access packets received by other nodes
	topo = new cTopology("topo");
	topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());


	VirtualMobilityManager *MOModule = dynamic_cast<VirtualMobilityManager*>
	(topo->getNode(s)->getModule()->getSubmodule("MobilityManager"));

	double s_x,s_y,s_z;

	s_x=MOModule->getLocation().x;
	s_y=MOModule->getLocation().y;
	s_z=MOModule->getLocation().z;
	
	VirtualMobilityManager *MOModule2 = dynamic_cast<VirtualMobilityManager*>
	(topo->getNode(d)->getModule()->getSubmodule("MobilityManager"));

	double d_x,d_y,d_z;

	d_x=MOModule2->getLocation().x;
	d_y=MOModule2->getLocation().y;
	d_z=MOModule2->getLocation().z;

	
	return sqrt((s_x-d_x)*(s_x-d_x)+(s_y-d_y)*(s_y-d_y)+(s_z-d_z)*(s_z-d_z));

}


vector<double>plnDir calPlnDir(int id){
	//adjust the n in pl formal
	cTopology *topo;	// temp variable to access packets received by other nodes
	topo = new cTopology("topo");
	topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());

	VirtualMobilityManager *MOModule = dynamic_cast<VirtualMobilityManager*>
	(topo->getNode(id)->getModule()->getSubmodule("MobilityManager"));

	double newx=MOModule->getLocation().x;
	double newy=MOModule->getLocation().y;
	double newz=MOModule->getLocation().z;
	//calculate the distance D across the body
	double distance_z=1-newz;
	double across_body_z=0;
	if(newz<0.8){
		across_body_z=0.2;
	}
	else{
		across_body_z=distance_z;
	}

	double rate=across_body_z/distance_z;

	double distance=sqrt((newx-1)*(newx-1)+(newy-1)*(newy-1)+(newz-1)*(newz-1));
	double across_body_distance=distance*rate;
	
	return {distance,across_body_distance};
}

double calPld(int id){
	//adjust the n in pl formal
	auto nums = calPlnDir(id);

	double distance=nums[0];
	double across_body_distance=nums[1];
	
	double newExponent=2.4;
	if(across_body_distance<0.4){
		newExponent=4;
	}else if(across_body_distance<0.6){
		newExponent=5;
	}else{
		newExponent=6;
	}

	double sig=4;
	if(id==1){
		sig=5.6;
	}else if(id==2||id==4||id==6){
		sig=3.38;
	}else if(id==3||id==7){
		sig=12.6;
	}else{
		sig=6.67;
	}

	return 55 + 10.0 * newExponent * log10(distance)+normal(0,sig);
}

void calN(int id){
	ofstream outfilee("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/n.txt", ios::app);

	auto nums = calPlnDir(id);

	double distance=nums[0];
	double across_body_distance=nums[1];

	double rate=across_body_distance/distance;
	double res=rate*7.4+(1-rate)*2;
	outfilee<<"node:"<<id<<"-"<<res<<endl;

}

void sortPld(vector<vector<double>> &pldSort){
					ofstream outfilee("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/pld.txt", ios::app);
					//outfilee<<"show pld:"<<endl;

	for(int id=1;id<8;id++){
		double pld=calPld(id);
		pldSort.push_back({pld,double(id)});
		outfilee<<pld<<' ';
		calN(id);
	}
	outfilee.close();
	sort(pldSort.begin(),pldSort.end());
}


vector<double> calDis(){
	double maxDis=0;
	vector<double> disArr;
	for(int i=1;i<8;i++){
		int tempDis=cal_dir(0,i);
		maxDis=tempDis>maxDis?tempDis:maxDis;
		disArr.push_back(tempDis);
	}
	for(int i=0;i<disArr.size();i++){
		disArr[i]/=maxDis;
	}
	return disArr;
}

void sortEcol(vector<vector<double>> &ecolSort){
	auto disArr=calDis();
	for(int id=1;id<8;id++){
		cTopology *topo;	// temp variable to access packets received by other nodes
		topo = new cTopology("topo");
		topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());
		ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
		(topo->getNode(id)->getModule()->getSubmodule("Application"));

		int energy_index=appModule->energy_index;
		auto energy_harvest=appModule->energy_harvest;
		double ecol=(energy_harvest[energy_index]/1000)*disArr[id-1];
		ecolSort.push_back({ecol,double(id)});
	}
	sort(ecolSort.begin(),ecolSort.end());
}


void ThroughputTest::startup()
{
	entropyData.resize(0);

	packet_rate = par("packet_rate");
	recipientAddress = par("nextRecipient").stringValue();
	recipientId = atoi(recipientAddress.c_str());
	startupDelay = par("startupDelay");
	delayLimit = par("delayLimit");
	packet_spacing = packet_rate > 0 ? 1 / float (packet_rate) : -1;
	
	cTopology *topo;		topo = new cTopology("topo");
		topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());
	ResourceManager *resmodule =  dynamic_cast<ResourceManager*>
				(topo->getNode(self)->getModule()->getSubmodule("ResourceManager"));

	double i_energy=resmodule->remainingEnergy;
	double max_energy=resmodule->initialEnergy;
	send_rate=base_send_rate+rand()%3;
	send_spacing=1/float(send_rate);
	dataSN = 0;
	
	numNodes = getParentModule()->getParentModule()->par("numNodes");
	packetsSent.clear();
	packetsReceived.clear();
	bytesReceived.clear();
	if(self==0)
	{
		resmodule->remainingEnergy*=200;
		resmodule->initialEnergy*=200;
		setTimer(CHANGE_SLOTS,startupDelay);
		setTimer(SHOW_ENERGY,startupDelay);
	}
	else{
		resmodule->remainingEnergy*=20;
		resmodule->initialEnergy*=20;
	}
	if(self==5)
	{
		setTimer(ENERGY_COL,startupDelay);
		//setTimer(SHOW_ENERGY,startupDelay);
		setTimer(WRITE_DATA,2);
		//setTimer(TEST_SLOTS,10);
		setTimer(ENERGY_TRAN,startupDelay);

	}
	if (self==7){
		setTimer(ENERGY_COL,startupDelay);
		setTimer(WRITE_DATA,2);
		setTimer(ENERGY_TRAN,startupDelay);
	}
	if (packet_spacing > 0 && recipientAddress.compare(to_string(self)) != 0)
	{
		setTimer(ADD_PACKET,packet_spacing+startupDelay);
		setTimer(SEND_PACKET, send_spacing + startupDelay);
	}
	else
		trace() << "Not sending packets";
	declareOutput("Packets received per node");
}

void ThroughputTest::fromNetworkLayer(ApplicationPacket * rcvPacket,
		const char *source, double rssi, double lqi)
{
	cTopology *topo;		
	topo = new cTopology("topo");
	topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());
	int sequenceNumber = rcvPacket->getSequenceNumber();
	int sourceId = atoi(source);

	// This node is the final recipient for the packet
	if (recipientAddress.compare(to_string(self)) == 0) {
		if (delayLimit == 0 || (simTime() - rcvPacket->getCreationTime()) <= delayLimit) { 
			trace() << "Received packet #" << sequenceNumber << " from node " << source;
			collectOutput("Packets received per node", sourceId);
			packetsReceived[sourceId]++;
			bytesReceived[sourceId] += rcvPacket->getByteLength();
		} else {
			trace() << "Packet #" << sequenceNumber << " from node " << source <<
				" exceeded delay limit of " << delayLimit << "s";
		}
	// Packet has to be forwarded to the next hop recipient
	} else {
		ApplicationPacket* fwdPacket = rcvPacket->dup();
		// Reset the size of the packet, otherwise the app overhead will keep adding on
		fwdPacket->setByteLength(0);
		toNetworkLayer(fwdPacket, recipientAddress.c_str());
	}
	//todo 异步处理packet_count--，利用timecallback即可
		ResourceManager *resmodule =  dynamic_cast<ResourceManager*>
				(topo->getNode(self)->getModule()->getSubmodule("ResourceManager"));
		double i_energy=resmodule->remainingEnergy;
		double max_energy=resmodule->initialEnergy;
		resmodule->consumeEnergy(13.36);
					trace()<<"node:"<<self<<"tx consume energy 6"<<endl;

		send_spacing=1/float(send_rate);
		setTimer(SINK_SEND,send_spacing);
		trace()<<"usleep："<<send_spacing<<"send_spacing"<<send_spacing<<endl;

}

void ThroughputTest::timerFiredCallback(int index)
{
	cTopology *topo;		topo = new cTopology("topo");
	topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());
	ResourceManager *resmodule =  dynamic_cast<ResourceManager*>
				(topo->getNode(self)->getModule()->getSubmodule("ResourceManager"));
	switch (index) {
		case SEND_PACKET:{
			sendCount++;
			if(self==0){
				packet_count=0;
				break;
			}
			double i_energy=resmodule->remainingEnergy;
			double max_energy=resmodule->initialEnergy;
			send_spacing=1/float(send_rate);

			double efs_mp=0;
			if(self!=0){
				double dir=cal_dir(0,self);
				if(dir<dmp){
					efs_mp=dir*dir*E_FS*800/1000;
				}else{
					efs_mp=dir*dir*dir*dir*E_MP*800/1000;
				}
			}
			resmodule->consumeEnergy(28.88+efs_mp);
			if(self==2)
				node_two_consume+=2;
			trace()<<"node:"<<self<<"tx consume energy:"<<28.88<<endl;
			trace() << "Sending packet #" << dataSN;
			toNetworkLayer(createGenericDataPacket(0, dataSN), recipientAddress.c_str());
			packetsSent[recipientId]++;
			dataSN++;
			//packet_count--;
		
			
			setTimer(SEND_PACKET, send_spacing);
			
				
			break;
		}
		// case ADD_PACKET:{
		// 	if(packet_count<packet_cap)
		// 		packet_count++;
		// 	setTimer(ADD_PACKET,packet_spacing);
		// 	break;
		// }
		case ENERGY_COL:{
			double temp;
			if(self==5){
				temp=energy_harvest[energy_index]/1000;
				energy_index++;
			}else{
				temp=energy_harvest[energy_index2]/1000;
				energy_index2++;
			}

			if(energy_index==10)
				energy_index=0;
			if(energy_index2==10)
				energy_index2=0;

			if(resmodule->remainingEnergy<=resmodule->initialEnergy-temp)
				resmodule->remainingEnergy+=temp;
			trace()<<"node:"<<self<<"col energy :"<<temp<<"remain: "<<resmodule->remainingEnergy<<endl;
			setTimer(ENERGY_COL,1);
			break;
		}
		case ENERGY_TRAN:{
			double i_energy=resmodule->remainingEnergy;
			double max_energy=resmodule->initialEnergy;
			double x,y,z;
			z=0;
			x=i_energy/max_energy;

			ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
					(topo->getNode(self)->getModule()->getSubmodule("Application"));
			y=float(appModule->packet_count)/appModule->packet_cap;
			for(int j=0;j<8;j++)
			{
				if(j==self)
					continue; 
				ThroughputTest *appModules = dynamic_cast<ThroughputTest*>
					(topo->getNode(j)->getModule()->getSubmodule("Application"));
				z+=(float(appModules->packet_count)/appModules->packet_cap);
			}
			z/=7;
			
			//double quanzhong=x+(1-y)+(1-z);
			
			double quanzhong;
			if(self==5){
				quanzhong=w[0]*x+(w[1]+w[2])*((1-y)/(1-z));
				max_energy_tx=energy_harvest[energy_index]/1000;

			}else{
				quanzhong=w2[0]*x+(w2[1]+w2[2])*((1-y)/(1-z));
				max_energy_tx=energy_harvest[energy_index2]/1000;

			}
			//trace()<<"send - Wwww:"<<w[0]<<' '<<w[1]<<' '<<w[2]<<endl;

			//double quanzhong=(x-y-z+2)/3;
			//double quanzhong=1;
			double energy_tran=max_energy_tx*quanzhong;
			resmodule->consumeEnergy(energy_tran*7/8);
			//node_two_consume+=energy_tran;
			trace()<<"node:"<<self<<"shuchu energy :"<<energy_tran<<" weight:"<<quanzhong<<"xyz"<<x<<' '<<y<<' '<<z<<' '<<endl;
			
			for(int j=0;j<8;j++)
			{
				if(j== 5||j==7)
					continue; 
				ResourceManager *rsmd =  dynamic_cast<ResourceManager*>
					(topo->getNode(j)->getModule()->getSubmodule("ResourceManager"));
				rsmd->remainingEnergy+=(energy_tran/6);
				trace()<<"node:"<<j<<"shuru energy :"<<energy_tran/6<<"from node"<<self<<"  "<<energy_tran<<" weight:"<<quanzhong<<"xyz"<<x<<' '<<y<<' '<<z<<' '<<endl;

			}
			
				
			
			setTimer(ENERGY_TRAN,1);
			break;

		}
		case SHOW_ENERGY:{
			for(int x=0;x<8;x++)
			{
				ResourceManager *resmodule =  dynamic_cast<ResourceManager*>
				(topo->getNode(x)->getModule()->getSubmodule("ResourceManager"));
				ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
					(topo->getNode(x)->getModule()->getSubmodule("Application"));
				trace()<<"node "<<x<<":"<<resmodule->remainingEnergy<<" count:"<<appModule->packet_count<<" all send:"<<appModule->sendCount<<endl;
			}
			//trace()<<"node5 consumeeee："<<node_two_consume<<endl;
			setTimer(SHOW_ENERGY,10);
			break;
		}
		case SINK_SEND:{
			// if(packet_count>0)
			// 	packet_count--;
			// break;
		}
		case WRITE_DATA:{
			
			if(self==5){
				double i_energy=resmodule->remainingEnergy;
				double max_energy=resmodule->initialEnergy;
				double x,y,z;
				z=0;
				x=float(i_energy)/max_energy;

				ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
						(topo->getNode(5)->getModule()->getSubmodule("Application"));
				y=float(appModule->packet_count)/appModule->packet_cap;
				for(int j=0;j<8;j++)
				{
					if(j==5)
						continue;
					ThroughputTest *appModules = dynamic_cast<ThroughputTest*>
						(topo->getNode(j)->getModule()->getSubmodule("Application"));
					z+=(float(appModules->packet_count)/appModules->packet_cap);
				}
				z/=7;
				entropyData.push_back({x+0.0001,y+0.0001,z+0.0001});
				if(entropyData.size()==100){
					//update weight
					w=entropyMethod(entropyData);
					
					for(int j=0;j<8;j++){
						if(j==5)
							continue;
						ThroughputTest *appModules = dynamic_cast<ThroughputTest*>
							(topo->getNode(j)->getModule()->getSubmodule("Application"));
						//trace()<<"gaibianqian:"<<appModules->w[0]<<endl;
						appModules->w=w;
						//trace()<<"gaibianhou:"<<appModules->w[0]<<endl;

					}
					trace()<<"change w1!:"<<w[0]<<' '<<w[1]<<' '<<w[2]<<endl;

					//resize datas
					entropyData.resize(0);
				}
			}
			if (self ==7){
				double i_energy=resmodule->remainingEnergy;
				double max_energy=resmodule->initialEnergy;
				double x,y,z;
				z=0;
				x=float(i_energy)/max_energy;

				ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
						(topo->getNode(7)->getModule()->getSubmodule("Application"));
				y=float(appModule->packet_count)/appModule->packet_cap;
				for(int j=0;j<8;j++)
				{
					if(j==7)
						continue;
					ThroughputTest *appModules = dynamic_cast<ThroughputTest*>
						(topo->getNode(j)->getModule()->getSubmodule("Application"));
					z+=(float(appModules->packet_count)/appModules->packet_cap);
				}
				z/=7;
				entropyData2.push_back({x+0.0001,y+0.0001,z+0.0001});
				if(entropyData2.size()==100){
					//update weight
					w2=entropyMethod(entropyData2);
					
					for(int j=0;j<8;j++){
						if(j==7)
							continue;
						ThroughputTest *appModules = dynamic_cast<ThroughputTest*>
							(topo->getNode(j)->getModule()->getSubmodule("Application"));
						//trace()<<"gaibianqian:"<<appModules->w[0]<<endl;
						appModules->w2=w2;
						//trace()<<"gaibianhou:"<<appModules->w[0]<<endl;

					}
					trace()<<"change Wwwww!:"<<w2[0]<<' '<<w2[1]<<' '<<w2[2]<<endl;

				//resize datas
				entropyData2.resize(0);
				}
			}
			
			setTimer(WRITE_DATA,2);
			break;
		}
		case CHANGE_SLOTS:{
			int packageCount=0;
			// calculate the number of packages in this round
			send_rate=base_send_rate+rand()%3;
			packageCount+=send_rate;
			for(int id=1;id<numNodes;id++){
				ThroughputTest *appmod = dynamic_cast<ThroughputTest*>
					(topo->getNode(id)->getModule()->getSubmodule("Application"));
				appmod->send_rate=base_send_rate+rand()%3;
				packageCount+=appmod->send_rate;
				trace()<<"node:"<<id<<"send_rate:"<<appmod->send_rate<<endl;
			}
			//trace()<<"this round all send_rate:"<<packageCount<<endl;
			vector <int> node_slots(numNodes,2);
				// ofstream outfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/gts", ios::trunc);
				// outfile << "2 2 2 2 2 2 2 2";
				
			//auction
			//packageCount=0;
			if(packageCount>packageCapacity){
				//this round dont neet auction
				vector<double> bid(numNodes,0);
				//evaluate
				//sort pld
				vector<vector<double>> pldSort;
				sortPld(pldSort);

				//sort energy harvest ability
				vector<vector<double>> ecolSort;
				sortEcol(ecolSort);

				map<int,int> mpPld;
				map<int,int> mpEcol;
				for(int index=0;index<pldSort.size();index++){
					mpPld.insert({int(pldSort[index][1]),index});
					mpEcol.insert({int(ecolSort[index][1]),index});
				}

				for(int id=0;id<numNodes;id++){
					ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
						(topo->getNode(id)->getModule()->getSubmodule("Application"));
					//appModule->money += (addMoney+rand()%20);
					int needTram=4+rand()%3;
					if(needTram==5){
						appModule->money+=5;
					}else if (needTram==6){
						appModule->money+=10;
					}

					appModule->money += addMoney;
					double base=1;

					ResourceManager *resmod =  dynamic_cast<ResourceManager*>
						(topo->getNode(id)->getModule()->getSubmodule("ResourceManager"));
					double energyRate = resmod->remainingEnergy / resmod->initialEnergy;

					if(energyRate>=0.6 && energyRate<0.8){
						base *= (energyRate+0.2);
					}else if(energyRate>=0.4 && energyRate<0.6){
						base *= (energyRate+0.1);
					}else if(energyRate>=0.2 && energyRate<0.4){
						base *= energyRate;
					}else if(energyRate < 0.2){
						base = 0.2;
					}

					double congestRate = double(appModule->packet_count)/double(appModule->packet_cap);
					if(congestRate>=0.05 && congestRate<0.1){
						base *= 0.5;
					}else if(congestRate>=0.1 && congestRate<0.15){
						base *= 0.7;
					}else if(congestRate>=0.15 && congestRate<0.2){
						base *= 0.9;
					}else if(congestRate >= 0.2){
						base *=1;
					}

					base *= (1-mpPld[id]*0.1);
					base *= (mpEcol[id]*0.1+0.3);

					if(base>1){
						base=1;
					}
					trace()<<"node-"<<id<<" base:"<<base<<endl;
					bid[id]=appModule->money * base;
				}
				
				for(int id=0;id<8;id++){
					trace()<<"bid: "<<bid[id]<<endl;
				}

				//distribute
				vector<vector<double>> node_bid;
				for(int id=0;id<numNodes;id++){
					node_bid.push_back({bid[id],id});
				}
				sort(node_bid.begin(),node_bid.end());
				//slots 3,3,2,2,2,2,1,1
				
				
				//pay
				vector<double> pay(numNodes,0);
				for(int k=0;k<numNodes;k++){
					int id=node_bid[k][1];

					if(k<2){
						node_slots[id]=4;
					}else if(k<6){
						node_slots[id]=8;
					}else{
						node_slots[id]=12;
					}
					
					
					ThroughputTest *appM = dynamic_cast<ThroughputTest*>
						(topo->getNode(id)->getModule()->getSubmodule("Application"));
					double needPay=0;
					if(k==0||k==1){
						needPay=node_bid[0][0];
					} else if (k>=2&&k<=5){
						needPay=node_bid[0][0]+node_bid[2][0]-node_bid[1][0];
					}else{
						needPay=node_bid[0][0]+node_bid[2][0]-node_bid[1][0]+node_bid[6][0]-node_bid[5][0];
					}
					pay[id]=needPay;

					appM->money -= needPay;
				}
				
				string slotStirng="";
				for(int k=0;k<numNodes;k++){
					if(k!=0){
						slotStirng+=" ";
					}
					slotStirng+=to_string(node_slots[k]);
				}
				ofstream outfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/gts", ios::trunc);
				outfile << slotStirng;
				outfile.close();

				for(int id=0;id<numNodes;id++){
					ThroughputTest *amd = dynamic_cast<ThroughputTest*>
						(topo->getNode(id)->getModule()->getSubmodule("Application"));
					if (node_slots[id]==4){
						int inbuffer= floor(0.35*amd->send_rate);
						amd->packet_count += inbuffer;
						amd->send_rate -= inbuffer;
					}else if(node_slots[id]==8){
						int inbuffer= floor(0.25*amd->send_rate);
						amd->packet_count += inbuffer;
						amd->send_rate -= inbuffer;
					}else{
					}
				
				}
				//test
				// for(int id=0;id<numNodes;id++){

				// 	ThroughputTest *aaa= dynamic_cast<ThroughputTest*>
				// 		(topo->getNode(id)->getModule()->getSubmodule("Application"));
				// 	trace()<<"node-"<<id<<" bid:"<<bid[id]<<"    pay:"<<pay[id]<<"   and remain money:"<<aaa->money<<endl;
				// }
			}else{
				//this round dont neet auction
				for(int id=0;id<numNodes;id++){
					ThroughputTest *app = dynamic_cast<ThroughputTest*>
						(topo->getNode(id)->getModule()->getSubmodule("Application"));
					
					int extra_send=floor(0.1*app->packet_count);
					// dont send too much, in order to decrease loss rate
					if(extra_send>15){
						extra_send=15;
					}
					app->packet_count -= extra_send;
					app->send_rate += extra_send;
				}
				ofstream outfile("/home/wangyiran/Castalia-master/Castalia/Simulations/BANtest/gts", ios::trunc);
				outfile << "8 8 8 8 8 8 8 8";
				outfile.close();
			}

			//output test
			trace()<<"this round slots:"<<endl;
			for(int id=0;id<numNodes;id++){
					ThroughputTest *am = dynamic_cast<ThroughputTest*>
						(topo->getNode(id)->getModule()->getSubmodule("Application"));
				trace()<<"node"<<id<<"send-rate"<<am->send_rate<<" slots:"<<node_slots[id]<<" buffer:"<<am->packet_count<<endl;
			
			}
			setTimer(CHANGE_SLOTS,1);
		}
		/*
		case TEST_SLOTS:{
			for(int id=0;id<5;id++){
				
				VirtualMac *test = dynamic_cast<VirtualMac*>
				(topo->getNode(id)->getModule()->getSubmodule("Communication")->getSubmodule("Mac"));
				trace()<<"node:"<<id<<" show gts list : "<<endl;
				for(auto k:test->GTSlist){
					trace()<<"slot owner:"<<k->owner<<" start:"<<k->start<<" len:"<<k->length<<endl;
				}
				test->requestGTS=2;
				for(auto k:test->GTSlist){
					trace()<<"slot owner:"<<k->owner<<" start:"<<k->start<<" len:"<<k->length<<endl;
				}
			}
						setTimer(TEST_SLOTS,10);


		}*/
	}

}


// This method processes a received carrier sense interupt. Used only for demo purposes
// in some simulations. Feel free to comment out the trace command.
void ThroughputTest::handleRadioControlMessage(RadioControlMessage *radioMsg)
{
	switch (radioMsg->getRadioControlMessageKind()) {
		case CARRIER_SENSE_INTERRUPT:
			trace() << "CS Interrupt received! current RSSI value is: " << radioModule->readRSSI();
                        break;
	}
}

void ThroughputTest::finishSpecific() {
	declareOutput("Packets reception rate");
	declareOutput("Packets loss rate");

	cTopology *topo;	// temp variable to access packets received by other nodes
	topo = new cTopology("topo");
	topo->extractByNedTypeName(cStringTokenizer("node.Node").asVector());

	long bytesDelivered = 0;
	for (int i = 0; i < numNodes; i++) {
		ThroughputTest *appModule = dynamic_cast<ThroughputTest*>
			(topo->getNode(i)->getModule()->getSubmodule("Application"));
		if (appModule) {
			int packetsSent = appModule->getPacketsSent(self);
			if (packetsSent > 0) { // this node sent us some packets
				float rate = (float)packetsReceived[i]/packetsSent;
				collectOutput("Packets reception rate", i, "total", rate);
				collectOutput("Packets loss rate", i, "total", 1-rate);
			}

			bytesDelivered += appModule->getBytesReceived(self);
		}
	}
	delete(topo);

	if (packet_rate > 0 && bytesDelivered > 0) {
		double energy = (resMgrModule->getSpentEnergy() * 1000000000)/(bytesDelivered * 8);	//in nanojoules/bit
		declareOutput("Energy nJ/bit");
		collectOutput("Energy nJ/bit","",energy);
	}
}



