import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import java.util.Scanner;
import java.util.concurrent.*;
import java.util.List;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import java.io.File;
import java.io.FileWriter;
import org.jlab.groot.ui.TGCanvas;

int num_elastic_events = 0; 
int countermax = 10000;
int num_files_max = (int)Double.parseDouble(args[1]);
double en = Double.parseDouble(args[0]);
double enmax = en+0.5; //GeV
double thetamax = 40;  //degrees
double phimax = 180;   //degrees
double vzmax = 50;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

File f = new File("output.txt");
f.delete();
FileWriter write = new FileWriter("output.txt",true);

System.out.println("Energy is " + en + " Max number of files is " + num_files_max);
HipoDataSource reader = new HipoDataSource();

H1F momentum = new H1F("momentum", "momentum", 500, 0, enmax);
momentum.setTitleX("momentum");

H1F W_hist = new H1F("W", "W", 500, 0, wmax);
W_hist.setTitleX("W");

H1F Q2_hist = new H1F("Q2", "Q2", 500, 0, enmax);
Q2_hist.setTitleX("Q2");

H2F W_vs_Q2 = new H2F("W_vs_Q2", "W_vs_Q2", 100, 0.0, enmax, 100, 0.0, wmax);
W_vs_Q2.setTitleX("Q2");
W_vs_Q2.setTitleY("W");

H2F E_vs_Theta = new H2F("E_vs_Theta", "E_vs_Theta", 100, 5, thetamax, 100, 0, enmax);
E_vs_Theta.setTitleX("Theta");
E_vs_Theta.setTitleY("E'");

H2F z_vs_Theta = new H2F("z_vs_Theta", "z_vs_Theta", 100, 5, thetamax, 100, -vzmax, vzmax);
z_vs_Theta.setTitleX("Theta");
z_vs_Theta.setTitleY("z vertex");

H2F Phi_vs_Theta = new H2F("Phi_vs_Theta", "Phi_vs_Theta", 100, 5, thetamax, 100, -phimax, phimax);
Phi_vs_Theta.setTitleX("Theta");
Phi_vs_Theta.setTitleY("Phi");

H2F Phi_vs_W = new H2F("Phi_vs_W", "Phi_vs_W", 100, 0, wmax, 100, -phimax, phimax);
Phi_vs_W.setTitleX("W");
Phi_vs_W.setTitleY("Phi");

H1F PhiResPE = new H1F("PhiResPE","PhiResPE", 500, 0, 360);
PhiResPE.setTitleX("Phi p - Phi e");

H1F EResPE = new H1F("EResPE","EResPE", 500, -10, 10);
EResPE.setTitleX("E p + E e - M p - Energy");

H1F VResPE = new H1F("VResPE","VResPE", 500, -50, 50);
VResPE.setTitleX("Vz P - Vz E");

H2F PhiPvsPhiE = new H2F("PhiPvsPhiE","PhiPvsPhiE",100, -180, 180, 100, -180, 180);
PhiPvsPhiE.setTitleX("Phi Electron");
PhiPvsPhiE.setTitleY("Phi Proton");

H2F EnPvsEnE = new H2F("EnPvsEnE","EnPvsEnE",100, 0 , enmax, 100, 0.9, 2);
EnPvsEnE.setTitleX("Energy Electron");
EnPvsEnE.setTitleY("Energy Proton");

H2F ZPvsZE = new H2F("ZPvsZE","ZPvsZE",100, -10 , 10, 100, -10, 10);
ZPvsZE.setTitleX("vertex Electron");
ZPvsZE.setTitleY("vertex Proton");

H1F momentum_el = new H1F("momentum_el", "momentum_el", 500, 0, enmax);
momentum_el.setTitleX("momentum");

H1F W_hist_el = new H1F("W_el", "W_el", 500, 0, 1.5);
W_hist_el.setTitleX("W");

H1F Q2_hist_el = new H1F("Q2_el", "Q2_el", 500, 0, enmax);
Q2_hist_el.setTitleX("Q2");

H2F W_vs_Q2_el = new H2F("W_vs_Q2_el", "W_vs_Q2_el", 100, 0.0, enmax, 100, 0.0, 1.5);
W_vs_Q2_el.setTitleX("Q2");
W_vs_Q2_el.setTitleY("W");

H2F E_vs_Theta_el = new H2F("E_vs_Theta_el", "E_vs_Theta_el", 100, 5, thetamax, 100, 0, enmax);
E_vs_Theta_el.setTitleX("Theta");
E_vs_Theta_el.setTitleY("E'");

H2F z_vs_Theta_el = new H2F("z_vs_Theta_el", "z_vs_Theta_el", 100, 5, thetamax, 100, -vzmax, vzmax);
z_vs_Theta_el.setTitleX("Theta");
z_vs_Theta_el.setTitleY("z vertex");

H2F Phi_vs_Theta_el = new H2F("Phi_vs_Theta_el", "Phi_vs_Theta_el", 100, 5, thetamax, 100, -phimax, phimax);
Phi_vs_Theta_el.setTitleX("Theta");
Phi_vs_Theta_el.setTitleY("Phi");

H2F Phi_vs_W_el = new H2F("Phi_vs_W_el", "Phi_vs_W_el", 100, 0, wmax, 100, -phimax, phimax);
Phi_vs_W_el.setTitleX("W");
Phi_vs_W_el.setTitleY("Phi");

H1F PhiResPE_el = new H1F("PhiResPE_el","PhiResPE_el", 500, 0, 360);
PhiResPE_el.setTitleX("Phi p - Phi e");

H1F EResPE_el = new H1F("EResPE_el","EResPE_el", 500, -10, 10);
EResPE_el.setTitleX("E p + E e - M p - Energy");

H1F VResPE_el = new H1F("VResPE_el","VResPE_el", 500, -50, 50);
VResPE_el.setTitleX("Vz P - Vz E");

H2F PhiPvsPhiE_el = new H2F("PhiPvsPhiE_el","PhiPvsPhiE_el",100, -180, 180, 100, -180, 180);
PhiPvsPhiE_el.setTitleX("Phi Electron");
PhiPvsPhiE_el.setTitleY("Phi Proton");

H2F EnPvsEnE_el = new H2F("EnPvsEnE_el","EnPvsEnE_el",100, 0 , enmax, 100, 0.9, 2);
EnPvsEnE_el.setTitleX("Energy Electron");
EnPvsEnE_el.setTitleY("Energy Proton");

H2F ZPvsZE_el = new H2F("ZPvsZE_el","ZPvsZE_el",100, -10 , 10, 100, -10, 10);
ZPvsZE_el.setTitleX("vertex Electron");
ZPvsZE_el.setTitleY("vertex Proton");

H2F Cal_y_vs_x_precut = new H2F("Cal_y_vs_x_precut", "Cal_y_vs_x_precut", 100, -450,450, 100, -450,450);
Cal_y_vs_x_precut.setTitleX("X (cm)");
Cal_y_vs_x_precut.setTitleY("Y (cm)");

H1F Cal_lu = new H1F("Cal_lu", "Cal_lu", 500, 0, 450); 
H1F Cal_lv = new H1F("Cal_lv", "Cal_lv", 500, 0, 450); 
H1F Cal_lw = new H1F("Cal_lw", "Cal_lw", 500, 0, 450); 

H2F Cal_y_vs_x = new H2F("Cal_y_vs_x", "Cal_y_vs_x", 100, -450,450, 100, -450, 450);
Cal_y_vs_x.setTitleX("X (cm)");
Cal_y_vs_x.setTitleY("Y (cm)");

H2F DCXvsHTCCX = new H2F("DCXvsHTCCX", "DCXvsHTCCX", 100, -450,450,100,-50,50);
DCXvsHTCCX.setTitleX("DC X");
DCXvsHTCCX.setTitleY("delta X");

H2F DCYvsHTCCY = new H2F("DCYvsHTCCY", "DCYvsHTCCY", 100, -450,450,100,-50,50);
DCYvsHTCCY.setTitleX("DC Y");
DCYvsHTCCY.setTitleY("delta Y");


H2F DCPhivsHTCCPhi = new H2F("DCPhivsHTCCPhi", "DCPhivsHTCCPhi", 100, -180,180,100,-10,10);
DCPhivsHTCCPhi.setTitleX("DC Phi");
DCPhivsHTCCPhi.setTitleY("delta Phi");

H2F ZPvsZEP95 = new H2F("ZPvsZEP95","ZPvsZEP95",100, -10 , 10, 100, -10, 10);
ZPvsZEP95.setTitleX("vertex Electron");
ZPvsZEP95.setTitleY("vertex Proton");

H1F PhiResPEP95 = new H1F("PhiResPEP95","PhiResPEP95", 500, 0, 360);
PhiResPEP95.setTitleX("Phi p - Phi e");

H1F EResPEP95 = new H1F("EResPEP95","EResPEP95", 500, -10, 10);
EResPEP95.setTitleX("E p + E e - M p - Energy");

H1F VResPEP95 = new H1F("VResPEP95","VResPEP95", 500, -50, 50);
VResPEP95.setTitleX("Vz P - Vz E");

H2F nphevsp = new H2F("nphevsp","hphevsp", 100, 0, 10, 100, 0, 50);
nphevsp.setTitleX("momentum electron");
nphevsp.setTitleY("n of photoelectrons");

H2F nphevstheta = new H2F("nphevstheta","nphevstheta",100,5,40,100,0,50);
nphevstheta.setTitleX("theta electron");
nphevstheta.setTitleY("n of photoelectrons");

H2F Cal_energy_p_vs_p = new H2F("cal_energy_p_vs_p","cal_energy_p_vs_p",100,0,10,100,0,0.5);
Cal_energy_p_vs_p.setTitleX("momentum electron");
Cal_energy_p_vs_p.setTitleY("Cal energy/p");

H2F ECInvsOut = new H2F("ECInvsOut","ECInvsOut",100,0,1,100,0,1);
ECInvsOut.setTitleX("EC Out");
ECInvsOut.setTitleY("EC In");

H2F ECvsPC = new H2F("ECvsPC","ECvsPC",100,0,1,100,0,1);
ECvsPC.setTitleX("PCal");
ECvsPC.setTitleY("ECal In+Out");


H1F p_hw_mom = new H1F("p_hw_mom","p_hw_mom",500,0,10);
p_hw_mom.setTitleX("Proton Momentum (W > 1.5)");

H1F p_hw_theta = new H1F("p_hw_theta","p_hw_theta",500,34,180);
p_hw_theta.setTitleX("Proton Theta (W > 1.5)");

H1F p_hw_phi = new H1F("p_hw_phi","p_hw_phi",500,-180,180);
p_hw_phi.setTitleX("Proton Phi (W > 1.5)");

H1F p_hw_momP95 = new H1F("p_hw_momP95","p_hw_momP95",500,0,10);
p_hw_momP95.setTitleX("Proton Momentum (W > 1.5)");

H1F p_hw_thetaP95 = new H1F("p_hw_thetaP95","p_hw_thetaP95",500,34,180);
p_hw_thetaP95.setTitleX("Proton Theta (W > 1.5)");

H1F p_hw_phiP95 = new H1F("p_hw_phiP95","p_hw_phiP95",500,-180,180);
p_hw_phiP95.setTitleX("Proton Phi (W > 1.5)");

H1F ECalc = new H1F("ECalc","ECalc",500,0,11);
ECalc.setTitleX("Energy Calculated using theta p and e"); 
	
H2F ThetaPvsThetaE = new H2F("ThetaPvsThetaE","ThetaPvsThetaE",100,0,20,100,0,180);
ThetaPvsThetaE.setTitleX("Theta Electron");
ThetaPvsThetaE.setTitleY("Theta Proton");

H1F ECalc_el = new H1F("ECalc_el","ECalc_el",500,0,11);
ECalc_el.setTitleX("Energy Calculated using theta p and e"); 
	
H2F ThetaPvsThetaE_el = new H2F("ThetaPvsThetaE_el","ThetaPvsThetaE_el",100,0,20,100,0,180);
ThetaPvsThetaE_el.setTitleX("Theta Electron");
ThetaPvsThetaE_el.setTitleY("Theta Proton");

H2F ProtonPvsTheta = new H2F("ProtonPvsTheta","ProtonPvsTheta",100,0,180,100,0,4);
ProtonPvsTheta.setTitleX("Theta");
ProtonPvsTheta.setTitleY("Proton Momentum");

H2F ProtonPvsTheta_el = new H2F("ProtonPvsTheta_el","ProtonPvsTheta_el",100,0,180,100,0,4);
ProtonPvsTheta_el.setTitleX("Theta");
ProtonPvsTheta_el.setTitleY("Proton Momentum");

H2F ProtonPvsTheta90 = new H2F("ProtonPvsTheta90","ProtonPvsTheta90",100,0,180,100,0,4);
ProtonPvsTheta90.setTitleX("Theta");
ProtonPvsTheta90.setTitleY("Proton Momentum");

H2F ProtonPvsThetaPQ = new H2F("ProtonPvsThetaPQ","ProtonPvsThetaPQ",100,0,180,100,0,4);
ProtonPvsThetaPQ.setTitleX("ThetaPQ");
ProtonPvsThetaPQ.setTitleY("Proton Momentum");

H1F ProtonP = new H1F("ProtonP", "ProtonP", 500, 0, 4);
ProtonP.setTitleX("Proton Momentum");

H1F ProtonP90 = new H1F("ProtonP90", "ProtonP90", 500, 0, 4);
ProtonP90.setTitleX("Proton Momentum for Theta > 90");

H2F ECalcvsTheta = new H2F("ECalcvsTheta","ECalcvsTheta",100,0,40,100,0,11);
ECalcvsTheta.setTitleX("Theta Electron");
ECalcvsTheta.setTitleY("Energy calculated using Theta proton and electron");

H2F ECalcvsTheta_el = new H2F("ECalcvsTheta_el","ECalcvsTheta_el",100,0,40,100,0,11);
ECalcvsTheta_el.setTitleX("Theta Electron");
ECalcvsTheta_el.setTitleY("Energy calculated using Theta proton and electron");


double e_mass = 0.000511;
double p_mass = 0.93827203;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);
int numfiles = args.length-2;
int count_success = 0; 
for(int i = 2; i < numfiles; i++){
if(i > num_files_max + 1){break;}
reader.open(args[i]);
double numevents = (double)reader.getSize();
int filenum = i-1;
double emax = 0;
phimax = 0;
thetamax = 0;
vzmax = 0;
double counter = 0;
byte sector = 0;
int cal_row = 0;
int dc_row = 0;
int htcc_row = 0;
float x_dc = 1000.0;
float y_dc = 1000.0;
float z_dc = 1000.0;
float x_htcc = 1000.0;
float y_htcc = 1000.0;
float z_htcc = 1000.0;
double phi_dc = 1000.0;
double phi_htcc = 1000.0;
float e_px = 0;
float e_py = 0;
float e_pz = 0; 
float e_beta = 0;
float e_vz = 0;
float e_phi = 0;
float e_mom = 0;
float e_theta = 0;

float p_px = 0;
float p_py = 0;
float p_pz = 0; 
float p_beta = 0;
float p_vz = 0;
float p_phi = 0;
float p_theta = 0; 
float p_mom = 0;
float pcal_energy = 0;
float ecin_energy = 0;
float ecout_energy = 0;
float tot_cal_energy = 0;
float lu = 0;
float lv = 0;
float lw = 0;
     
double[] ctof_info = new double[2];

ConcurrentHashMap<Integer,Float> cal_energy_map = new ConcurrentHashMap<Integer,Float>();

float x_cal = 0;
float y_cal = 0;
double W = 0; 
double Q2 = 0;


int htccrow = -1;
double nphe = 0;
double htcc_theta = 0;
boolean ec_flag = false;
boolean found_electron = false;
boolean found_proton = false; 
boolean dcxhtccx_flag = false;
boolean dcyhtccy_flag = false;
boolean dcphihtccphi_flag = false;
boolean dc_flag = false;
boolean htcc_flag = false;
boolean pid_flag = false;
boolean dt_flag = false; 
Vector3 e_vec_3 = new Vector3();
Vector3 p_vec_3 = new Vector3();
LorentzVector e_vec_prime = new LorentzVector();
LorentzVector p_vec_prime = new LorentzVector();
LorentzVector q_vec = new LorentzVector();
double theta_dc = 0;
double pos_dc = 0;
double theta_htcc = 0;
double pos_htcc = 0;
double path_length = 0; 
double dt = 0; 
double time_ctof = 0; 
while (reader.hasEvent()) {
        ec_flag = false;
	dc_flag = false; 
	pid_flag = false;
	htcc_flag = false; 
	dcxhtccx_flag = false;
	dcyhtccy_flag = false;
	dt_flag = false;
	found_electron = false;
	found_proton = false;
	dcphihtccphi_flag = false;
	DataEvent event = reader.getNextEvent();
	counter++;

	System.out.print("file " + filenum + "/" + numfiles + " "  + counter*100/numevents + "% complete.\r");
	if (event.hasBank("REC::Particle") && event.hasBank("REC::Cherenkov") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Traj")) {
		DataBank bank_rec = event.getBank("REC::Particle");
		DataBank bank_cal = event.getBank("REC::Calorimeter");
		DataBank bank_traj = event.getBank("REC::Traj");
		DataBank bank_htcc = event.getBank("REC::Cherenkov");

			
			//counter++;
			//if(counter > countermax){break;}
		for (int k = 0; k < bank_rec.rows(); k++) {
			int pid = bank_rec.getInt("pid", k);
			byte q = bank_rec.getByte("charge", k);
			float px = bank_rec.getFloat("px", k);
			float py = bank_rec.getFloat("py", k);								
			float pz = bank_rec.getFloat("pz", k);
			float beta = bank_rec.getFloat("beta", k);
			float mom = (float) Math.sqrt(px * px + py * py + pz * pz);
			double phi = Math.atan2((double) py,(double) px);
			double theta = Math.acos((double) pz/(double) mom);
			theta *= 180/Math.PI;
			phi *= 180/Math.PI;
			float vz = bank_rec.getFloat("vz", k);	     

			if(pid == 11 && !found_electron){
				
				e_px = px; 
				e_py = py; 
				e_pz = pz; 
				e_beta = beta;
				e_mom = mom; 
				e_phi = phi; 
				e_theta = theta; 
				e_vz = vz; 
				
				e_vec_3 = new Vector3(px, py, pz); //3 vector e'
				e_vec_prime = new LorentzVector(); //4 vector e'
				e_vec_prime.setVectM(e_vec_3, e_mass);

				if((e_vec_prime.e() > 0.1 * en) && theta > 5 && theta < 40){
					  found_electron = true;
					  q_vec = new LorentzVector(); //4 vector q
					  q_vec.copy(e_vec); //e - e'
					  q_vec.sub(e_vec_prime);
					  Q2 = -q_vec.mass2(); //-q^2
				
					  LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
					  w_vec.copy(p_vec); //p-q
					  w_vec.add(q_vec);
					  W = w_vec.mass(); 					 
					  htccrow = get_htcc_row(event,k);	
					  nphe = bank_htcc.getFloat("nphe",htccrow);
					  htcc_theta = bank_htcc.getFloat("theta",htccrow);		   
				          cal_energy_map = cal_cut_row(event,k);		

						if(!ec_flag && cal_energy_map.size() == 3){
								
							pcal_energy = 0;
							ecin_energy = 0;
							ecout_energy = 0;
							tot_cal_energy = 0;
								
							sector = bank_cal.getByte("sector",cal_row);
							x_cal = bank_cal.getFloat("x",cal_row);
							y_cal = bank_cal.getFloat("y",cal_row);
								
							Cal_y_vs_x_precut.fill(x_cal,y_cal);
							lu = cal_energy_map.get(1).get(1);
							lv = cal_energy_map.get(1).get(2);
							lw = cal_energy_map.get(1).get(3);
							pcal_energy = cal_energy_map.get(1).get(0);
							ecin_energy = cal_energy_map.get(4).get(0);
							ecout_energy = cal_energy_map.get(7).get(0);
 
							//System.out.println(ecin_energy + " " + ecout_energy + "      ");
							Cal_lu.fill(lu);
							Cal_lv.fill(lv);
							Cal_lw.fill(lw);				       
							if(ecin_energy > 0.06){// && lu > 60 && lu < 350 && lv < 370 && lw < 360){
								ec_flag = true;
								tot_cal_energy = pcal_energy + ecin_energy + ecout_energy;
							}	
						}
		
				} 
			}
				
			dc_row = dc_cut_row(event, k);
			if(dc_row != -1 && !dc_flag){
				x_dc = bank_traj.getFloat("x",dc_row);
				y_dc = bank_traj.getFloat("y",dc_row);
				z_dc = bank_traj.getFloat("z",dc_row);
				pos_dc = Math.sqrt(x_dc*x_dc + y_dc*y_dc + z_dc*z_dc);
				theta_dc = Math.acos((double) z_dc/ pos_dc);
				phi_dc = Math.atan2((double) y_dc,(double) x_dc);
				theta_dc *= 180/Math.PI;
				phi_dc *= 180/Math.PI;
				dc_flag = true;
			}
			
			htcc_row = htcc_cut_row(event,k);
			if(htcc_row != -1 && !htcc_flag){
				x_htcc = bank_traj.getFloat("x",htcc_row);
				y_htcc = bank_traj.getFloat("y",htcc_row);
				z_htcc = bank_traj.getFloat("z",htcc_row);
				pos_htcc = Math.sqrt(x_htcc*x_htcc + y_htcc*y_htcc + z_htcc*z_htcc);
				theta_htcc = Math.acos((double) z_htcc/ pos_htcc);
				phi_htcc = Math.atan2((double) y_htcc,(double) x_htcc);
				theta_htcc *= 180/Math.PI;
				phi_htcc *= 180/Math.PI;
				htcc_flag = true; 
			}

			if(pid == 2212 && !found_proton){
			       	found_proton = true;
				p_px = px; 
				p_py = py; 
				p_pz = pz; 
				p_beta = beta;
				p_mom = mom; 
				p_phi = phi; 
				p_theta = theta; 
				p_vz = vz; 
				
				p_vec_3 = new Vector3(px, py, pz); //3 vector e'
				p_vec_prime = new LorentzVector(); //4 vector e'
				p_vec_prime.setVectM(p_vec_3, p_mass);
				
				}

			}

		}
		if(found_electron && found_proton && nphe > 2){
		        ProtonP.fill(p_mom);
			ProtonPvsTheta.fill(p_theta,p_mom);
			ProtonPvsThetaPQ.fill(p_vec_3.theta(q_vec.vect()),p_mom);
		        ThetaPvsThetaE.fill(e_theta,p_theta);
			ECalc.fill(p_mass*(1/Math.tan(p_theta*Math.PI/180)/Math.tan(e_theta*Math.PI/360) - 1));
			if(p_theta > 90){
				   ProtonPvsTheta90.fill(p_theta,p_mom);
				   ProtonP90.fill(p_mom);
				   EnPvsEnE.fill(e_vec_prime.e(),p_vec_prime.e());
				   ZPvsZE.fill(e_vz,p_vz);
				   double eres = e_vec_prime.e()+p_vec_prime.e()-p_mass-en;
				   EResPE.fill(-eres);
				   VResPE.fill(e_vz-p_vz);
			           momentum.fill(e_mom);
				   W_hist.fill(W);
				   W_vs_Q2.fill(Q2,W);
				   Phi_vs_W.fill(W,e_phi);
				   Q2_hist.fill(Q2);
				   E_vs_Theta.fill(e_theta,e_vec_prime.e());
				   z_vs_Theta.fill(e_theta,e_vz);
				   Phi_vs_Theta.fill(e_theta,e_phi);
				   if(W > 1.5){
			     	   	p_hw_mom.fill(p_mom);
			     		p_hw_theta.fill(p_theta);
			     		p_hw_phi.fill(p_phi);
			     		if(p_theta > 95){
			     			   p_hw_momP95.fill(p_mom);
						   p_hw_thetaP95.fill(p_theta);
						   p_hw_phiP95.fill(p_phi);	
			     		}
			}
			}
		}		
		if(ec_flag && found_electron && found_proton && nphe > 2){
			PhiResPE.fill(Math.abs(e_phi-p_phi));
			PhiPvsPhiE.fill(e_phi,p_phi);
		}
		
//		System.out.println(Boolean.toString(ec_flag) + " " + Boolean.toString(found_electron) + " " + Boolean.toString(found_proton) + " " + Math.abs(e_phi-p_phi) + " " + nphe);		
		
		if(ec_flag && p_theta > 40 && found_electron && found_proton && nphe > 2){
			
			double eres = e_vec_prime.e()+p_vec_prime.e()-p_mass-en;
			//if(htccrow != -1){
			nphevsp.fill(e_mom,nphe);
			nphevstheta.fill(e_theta,nphe);
			//Cal_lu.fill(lu);
			//Cal_lv.fill(lv);
			//Cal_lw.fill(lw);
			Cal_y_vs_x.fill(x_cal,y_cal);

			Cal_energy_p_vs_p.fill(e_mom,tot_cal_energy/e_mom);
			ECInvsOut.fill(ecout_energy,ecin_energy);
			ECvsPC.fill(pcal_energy,ecin_energy+ecout_energy);

			if(dc_flag && htcc_flag){ 
				   DCXvsHTCCX.fill(x_dc,x_dc-x_htcc);
				   DCYvsHTCCY.fill(y_dc,y_dc-y_htcc);
				   DCPhivsHTCCPhi.fill(phi_dc,phi_dc-phi_htcc);
			}
			
						
			if(p_theta > 95){
			     PhiResPEP95.fill(Math.abs(e_phi-p_phi));
			     EResPEP95.fill(-eres);
			     VResPEP95.fill(e_vz-p_vz);
			     ZPvsZEP95.fill(e_vz,p_vz);
			}
			if(Math.abs(e_phi-p_phi) < 185 && Math.abs(e_phi-p_phi) > 175){
				EResPE_el.fill(-eres);
			}
			if(Math.abs(e_phi-p_phi) < 185 && Math.abs(e_phi-p_phi) > 175 && Math.abs(eres) <= 0.5){
				EnPvsEnE_el.fill(e_vec_prime.e(),p_vec_prime.e());
				
				momentum_el.fill(e_mom);
				W_hist_el.fill(W);
				W_vs_Q2_el.fill(Q2,W);
				Phi_vs_W_el.fill(W,e_phi);
				Q2_hist_el.fill(Q2);
				E_vs_Theta_el.fill(e_theta,e_vec_prime.e());
				ProtonPvsTheta_el.fill(p_theta,p_mom);
				if(W<=1){
					num_elastic_events++;
					PhiPvsPhiE_el.fill(e_phi,p_phi);
					PhiResPE_el.fill(Math.abs(e_phi - p_phi));
					ZPvsZE_el.fill(e_vz,p_vz);
					VResPE_el.fill(e_vz-p_vz);
					z_vs_Theta_el.fill(e_theta,e_vz);
					Phi_vs_Theta_el.fill(e_theta,e_phi);
					ThetaPvsThetaE_el.fill(e_theta,p_theta);
					ECalc_el.fill(p_mass*(1/Math.tan(p_theta*Math.PI/180)/Math.tan(e_theta*Math.PI/360) - 1));
				}
			} 
		}
		found_electron = false; 
		found_proton = false;
	}
}

boolean dc_cut(float X, float Y, int S)
{ 
	boolean result= false;
	if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
	&& (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
  
	return result;
}

int get_htcc_row(DataEvent event, int row){
       DataBank bank_htcc = event.getBank("REC::Cherenkov");
       int row_index = 0;
       for(int j = 0; j < bank_htcc.rows(); j++){
       	       row_index = bank_htcc.getInt("pindex",j);
	       if(row_index == row){
	       		    return j;	       
	       }
       }
       return -1;
}

ConcurrentHashMap<Integer,List<Float>> cal_cut_row(DataEvent event, int row){
	DataBank bank_cal = event.getBank("REC::Calorimeter");
	int row_index = 0;
	int cal_row_match = -1;
	int n_layers = 0;
	int layer = 0; 
	ConcurrentHashMap<Integer,List<Float>> cal_energy_map = new ConcurrentHashMap<Integer,List<Float>>();
	List<Float> l = new ArrayList<Float>();
	for(int j = 0; j < bank_cal.rows(); j++){
		row_index = bank_cal.getInt("pindex",j);
		if(row_index == row){
			l.add(0,bank_cal.getFloat("energy",j));
			l.add(1,bank_cal.getFloat("lu",j));
			l.add(2,bank_cal.getFloat("lv",j));
			l.add(3,bank_cal.getFloat("lw",j));
			layer = bank_cal.getInt("layer",j);
			cal_energy_map.put(layer,l);
		}
	}
	return cal_energy_map;
}

int dc_cut_row(DataEvent event, int row){
	DataBank bank_traj = event.getBank("REC::Traj");
	int row_index = 0;
	int det_id = 0;
	int cal_row_match = -1;
	float path_length = 0;
	for(int j = 0; j < bank_traj.rows(); j++){
		row_index = bank_traj.getInt("pindex",j);
		det_id = bank_traj.getInt("detId",j);
		path_length = bank_traj.getFloat("pathlength",j);
		if(row_index == row && path_length > 250 && path_length < 550){
			cal_row_match = j;
			break;
		}
	}
	return cal_row_match;
}
int htcc_cut_row(DataEvent event, int row){
	DataBank bank_traj = event.getBank("REC::Traj");
	int row_index = 0;
	int det_id = 0;
	int htcc_row_match = -1;
	float path_length = 0;
	for(int j = 0; j < bank_traj.rows(); j++){
		row_index = bank_traj.getInt("pindex",j);
		det_id = bank_traj.getInt("detId",j);
		//path_length = bank_traj.getFloat("pathlength",j);
		if(row_index == row && det_id == 0){
			htcc_row_match = j;
			break;
		}
	}
	return htcc_row_match;
}


//System.out.println(emax + " " + thetamax + " " + phimax + " " + vzmax);

/*std::cout << "\nNext Event?" << std::endl;                    
int key = std::cin.get();    
if (key == EOF || key == 'n' || key == 'N') return 0;                          
if (key != '\n') std::cin.ignore(numeric_limits<streamsize>::max(), '\n');     
std::cout << "OK" << std::endl;
*/

/*Scanner scanner = new Scanner(System.in);
System.out.println("press a key to continue");
scanner.nextLine();
System.out.println("OK");
*/

TGCanvas c1 = new TGCanvas("c_ep_coincidence","c_ep_coincidence",800,800);
c1.addCanvas("ProtonPvsTheta");
c1.draw(ProtonPvsTheta);

c1.addCanvas("ProtonPvsThetaPQ");
c1.draw(ProtonPvsThetaPQ);

c1.addCanvas("ProtonP");
c1.draw(ProtonP);

c1.addCanvas("ThetaPvsThetaE");
c1.draw(ThetaPvsThetaE);

c1.addCanvas("ECalc");
c1.draw(ECalc);


TGCanvas c2 = new TGCanvas("c_ep_coincidence_ptheta90","c_ep_coincidence_ptheta90",800,800);
c2.addCanvas("ProtonP90");
c2.draw(ProtonP90);

c2.addCanvas("EnPvsEnE");
c2.draw(EnPvsEnE);

c2.addCanvas("ZPvsZE");
c2.draw(ZPvsZE);

c2.addCanvas("EResPE");
c2.draw(EResPE);

c2.addCanvas("VResPE");
c2.draw(VResPE);

c2.addCanvas("momentum");
c2.draw(momentum);

c2.addCanvas("W_vs_Q2");
c2.draw(W_vs_Q2);

c2.addCanvas("Q2_hist");
c2.draw(Q2_hist);

c2.addCanvas("W_hist");
c2.draw(W_hist);

c2.addCanvas("E_vs_Theta");
c2.draw(E_vs_Theta);

c2.addCanvas("z_vs_Theta");
c2.draw(z_vs_Theta);

c2.addCanvas("Phi_vs_Theta");
c2.draw(Phi_vs_Theta);

c2.addCanvas("Phi_vs_W");
c2.draw(Phi_vs_W);

TGCanvas c3 = new TGCanvas("c_ep_coincidence_ptheta90_W1.5","c_ep_coincidence_ptheta90_W1.5",800,800);
c3.addCanvas("p_hw_mom");
c3.draw(p_hw_mom);

c3.addCanvas("p_hw_theta");
c3.draw(p_hw_theta);

c3.addCanvas("p_hw_phi");
c3.draw(p_hw_phi);

TGCanvas c4 = new TGCanvas("c_ep_coincidence_ptheta95_W1.5","c_ep_coincidence_ptheta95_W1.5",800,800);
c4.addCanvas("p_hw_momP95");
c4.draw(p_hw_momP95);

c4.addCanvas("p_hw_thetaP95");
c4.draw(p_hw_thetaP95);

c4.addCanvas("p_hw_phiP95");
c4.draw(p_hw_phiP95);

TGCanvas c5 = new TGCanvas("c_ec_epcoincidence","c_ec_epcoincidence",800,800);
c5.addCanvas("PhiPvsPhiE");
c5.draw(PhiPvsPhiE);

c5.addCanvas("PhiResPE");
c5.draw(PhiResPE);

TGCanvas c6 = new TGCanvas("c_ec_epcoincidence_ptheta40","c_ec_epcoincidence_ptheta40",800,800);
c6.addCanvas("nphevsp");
c6.draw(nphevsp);

c6.addCanvas("nphevstheta");
c6.draw(nphevstheta);

c6.addCanvas("Cal_energy_p_vs_p");
c6.draw(Cal_energy_p_vs_p);

c6.addCanvas("ECInvsOut");
c6.draw(ECInvsOut);

c6.addCanvas("ECvsPC");
c6.draw(ECvsPC);

c6.addCanvas("Cal_y_vs_x_precut");
c6.draw(Cal_y_vs_x_precut);

c6.addCanvas("Cal_y_vs_x");
c6.draw(Cal_y_vs_x);

c6.addCanvas("Cal_lu");
c6.draw(Cal_lu);

c6.addCanvas("Cal_lv");
c6.draw(Cal_lv);

c6.addCanvas("Cal_lw");
c6.draw(Cal_lw);

TGCanvas c7 = new TGCanvas("c_ec_epcoincidence_ptheta95","c_ec_epcoincidence_ptheta95",800,800);
c7.addCanvas("EResPEP95");
c7.draw(EResPEP95);

c7.addCanvas("VResPEP95");
c7.draw(VResPEP95);

c7.addCanvas("ZPvsZEP95");
c7.draw(ZPvsZEP95);

c7.addCanvas("PhiResPEP95");
c7.draw(PhiResPEP95);



TGCanvas canall = new TGCanvas("canall","canall", 800,800);


canall.addCanvas("DCXvsHTCCX");
canall.draw(DCXvsHTCCX);

canall.addCanvas("DCYvsHTCCY");
canall.draw(DCYvsHTCCY);

canall.addCanvas("DCPhivsHTCCPhi");
canall.draw(DCPhivsHTCCPhi);




TGCanvas c8 = new TGCanvas("c_ec_epcoincidence_ptheta40_phiresd5_eres<0.5","c_ec_epcoincidence_ptheta40_phiresd5_eres<0.5", 800,800);
c8.addCanvas("W_vs_Q2_el");
c8.draw(W_vs_Q2_el);

c8.addCanvas("momentum_el");
c8.draw(momentum_el);

c8.addCanvas("Q2_hist_el");
c8.draw(Q2_hist_el);

c8.addCanvas("W_hist_el");
c8.draw(W_hist_el);

c8.addCanvas("E_vs_Theta_el");
c8.draw(E_vs_Theta_el);

c8.addCanvas("Phi_vs_W_el");
c8.draw(Phi_vs_W_el);

c8.addCanvas("EnPvsEnE_el");
c8.draw(EnPvsEnE_el);

c8.addCanvas("PhiResPE_el");
c8.draw(PhiResPE_el);

c8.addCanvas("EResPE_el");
c8.draw(EResPE_el);



TGCanvas c9 = new TGCanvas("c_ec_epcoincidence_ptheta40_phiresd5_eres<0.5_W<1","c_ec_epcoincidence_ptheta40_phiresd5_eres<0.5_W<1", 800,800);

c9.addCanvas("PhiPvsPhiE_el");
c9.draw(PhiPvsPhiE_el);

c9.addCanvas("ZPvsZE_el");
c9.draw(ZPvsZE_el);

c9.addCanvas("VResPE_el");
c9.draw(VResPE_el);

c9.addCanvas("z_vs_Theta_el");
c9.draw(z_vs_Theta_el);

c9.addCanvas("Phi_vs_Theta_el");
c9.draw(Phi_vs_Theta_el);

c9.addCanvas("ThetaPvsThetaE_el");
c9.draw(ThetaPvsThetaE_el);

c9.addCanvas("ECalc_el");
c9.draw(ECalc_el);


System.out.println("Done!");
