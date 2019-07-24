import java.lang.reflect.Method;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.FunctionFactory;
import org.jlab.groot.math.UserParameter;
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
import org.jlab.groot.data.TDirectory
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.io.HipoWriter
TDirectory dir = new TDirectory();


double num_elastic_events = 0; 
int countermax = 10000;
boolean usecounter = true;
//boolean usecounter = true;
//boolean usepercentdisp = false;
boolean usepercentdisp = true;
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
//HipoDataSource reader = new HipoDataSource();
HipoReader reader = new HipoReader();
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

float[] bins = new float[16];//{0.18f,0.21f,0.25f,0.31f,0.36f,0.43f,0.51f,0.62f,0.73f,0.88f,1.05f,1.25f,1.49f,1.78f,2.12f,2.54f};
float bin_edge = 0.92;
float bin_edge_prev = 0; 
float bin_size = 0;
int index = 0; 
HashMap<Integer,Float> binmap = new HashMap<Integer,Float>();
HashMap<Integer,Integer> binmapcount = new HashMap<Integer,Integer>();
while(bin_edge < 15.64){
	       bin_edge_prev = bin_edge;
	       bins[index] = bin_edge;
	       binmap.put(index,bin_edge);
	       binmapcount.put(index,0);
	       bin_edge*=Math.pow((double)10,(double)1/13)	   
	       bin_size = bin_edge - bin_edge_prev;
	       index++;
}
for(int key : binmap.keySet()){
	System.out.println("debug " + binmap.get(key));
}
H1F Q2_hist_el = new H1F("Q2_el", 0.92, 15.63, bins);
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

/*H2F DCXvsHTCCX = new H2F("DCXvsHTCCX", "DCXvsHTCCX", 100, -450,450,100,-50,50);
DCXvsHTCCX.setTitleX("DC X");
DCXvsHTCCX.setTitleY("delta X");

H2F DCYvsHTCCY = new H2F("DCYvsHTCCY", "DCYvsHTCCY", 100, -450,450,100,-50,50);
DCYvsHTCCY.setTitleX("DC Y");
DCYvsHTCCY.setTitleY("delta Y");


H2F DCPhivsHTCCPhi = new H2F("DCPhivsHTCCPhi", "DCPhivsHTCCPhi", 100, -180,180,100,-10,10);
DCPhivsHTCCPhi.setTitleX("DC Phi");
DCPhivsHTCCPhi.setTitleY("delta Phi");
*/
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

H1F xb = new H1F("xb","xb",500,0,1);
xb.setTitleX("x");
H1F xbc = new H1F("xbc","xbc",500,0,1);
xbc.setTitleX("xc");
H2F xbcvsx = new H2F("xbcvsx","xbcvsx",100,0,1,100,0,1);
xbcvsx.setTitleX("x");
xbcvsx.setTitleY("y");

H1F[] W_pbins = new H1F[6];
H1F[] Q_pbins = new H1F[6];
H1F[] Pprot_pbins = new H1F[6];
H2F[] PprotvsThetaprot_pbins = new H2F[6];
H1F[] Wcorr_pbins = new H1F[6];
H1F[] Wcorr2_pbins = new H1F[6];
for(int i = 0; i < 6; i++){
	W_pbins[i] = new H1F("W_pbins_" + i,"W_pbins_" + i,500,0,wmax);
	W_pbins[i].setTitleX("W");
	Q_pbins[i] = new H1F("Q_pbins_" + i,"Q_pbins_" + i,500,0,enmax);
	Q_pbins[i].setTitleX("Q2");
	Pprot_pbins[i] = new H1F("Pprot_pbins_" + i,"Pprot_pbins_" + i,500,0,4);
	Pprot_pbins[i].setTitleX("Proton Momentum");
	PprotvsThetaprot_pbins[i] = new H2F("PprotvsThetaprot_pbins_" + i,"PprotvsThetaprot_pbins_" + i,100,89,181,100,0,4);
	PprotvsThetaprot_pbins[i].setTitleX("Proton Theta");
	PprotvsThetaprot_pbins[i].setTitleY("Proton Momentum");
	Wcorr_pbins[i] = new H1F("Wcorr_pbins_" + i,"Wcorr_pbins_" + i,500,0,wmax);
	Wcorr_pbins[i].setTitleX("W corrected (1.876)");
	Wcorr2_pbins[i] = new H1F("Wcorr2_pbins_" + i,"Wcorr2_pbins_" + i,500,0,wmax);
	Wcorr2_pbins[i].setTitleX("W corrected (1.86)");
}



double numevents = 0; 
double e_mass = 0.000511;
double p_mass = 0.93827203;
double n_mass = 0.93957;
double d_mass = 1.876;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);
int numfiles = args.length-2;
int count_success = 0; 
double counter = 0;
double banks_counter = 0; 
double ptheta_epcoin_nphe_count = 0;
double phi_count = 0;
double phi_eres_count = 0;
double w_count = 0;
for(int i = 2; i <= numfiles+1; i++){
	if(i > num_files_max + 1 && num_files_max != 0){break;}
	reader.open(args[i]);
	counter = 0;
	banks_counter = 0; 
	ec_ptheta_epcoin_nphe_count = 0;
	phi_eres_count = 0;
	w_count = 0;
	/*Method[] methods = reader.class.getDeclaredMethods();
	int nMethod = 1;
	System.out.println("1. List of all methods of Person class");
	for (Method method : methods) {
	System.out.printf("%d. %s", ++nMethod, method);
	System.out.println();
	}*/
	numevents = (double)reader.getEventCount();
	//numevents = 23;
	int filenum = i-1;
	double emax = 0;
	phimax = 0;
	thetamax = 0;
	vzmax = 0;

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
	double elec_e = 0;
	double elec_eprime = 0; 
	Vector3 q_vec_3 = new Vector3();
	Bank bank_rec = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
	Bank bank_cal = new Bank(reader.getSchemaFactory().getSchema("REC::Calorimeter"));
	Bank bank_traj = new Bank(reader.getSchemaFactory().getSchema("REC::Traj"));
	Bank bank_htcc = new Bank(reader.getSchemaFactory().getSchema("REC::Cherenkov"));
	Event event = new Event();


	while (reader.hasNext()) {
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
		//DataEvent event = reader.getNextEvent();
		reader.nextEvent(event);
		event.read(bank_rec);
		event.read(bank_cal);
		event.read(bank_traj);
		event.read(bank_htcc);
		counter++;

		if(usepercentdisp) System.out.print("file " + filenum + "/" + numfiles + " "  + counter*100/numevents + "% complete.\r");
		/****
		*Physics Calculations
		****/
		if(true){//if (event.hasBank("REC::Particle") && event.hasBank("REC::Cherenkov") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Traj")) {
			banks_counter++;
			//DataBank bank_rec = event.getBank("REC::Particle");

			//DataBank bank_cal = event.getBank("REC::Calorimeter");
			//DataBank bank_traj = event.getBank("REC::Traj");
			//DataBank bank_htcc = event.getBank("REC::Cherenkov");

			Map<Integer,List<Integer>> calMap = loadMapByIndex(bank_cal,"pindex");
			Map<Integer,List<Integer>> trajMap = loadMapByIndex(bank_traj,"pindex");
			Map<Integer,List<Integer>> htccMap = loadMapByIndex(bank_htcc,"pindex");



			//counter++;
			if(counter > countermax && usecounter){break;}
			for (int k = 0; k < bank_rec.getRows(); k++) {
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

				if(pid == 11 && k == 0 && !found_electron){

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
						q_vec_3 = q_vec.vect();
						Q2 = -q_vec.mass2(); //-q^2
						elec_e = en; 
						elec_eprime = e_vec_prime.e();
						LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
						w_vec.copy(p_vec); //p-q
						w_vec.add(q_vec);
						W = w_vec.mass(); 					 
						//htccrow = get_htcc_row(event,k);	
						for(int ihtcc : htccMap.get(k)){
							nphe = bank_htcc.getFloat("nphe",ihtcc);
							//htcc_theta = bank_htcc.getFloat("theta",ihtcc);		   
						}
						//System.out.println(nphe + " " + bank_htcc.getFloat("nphe",htccrow));
						//nphe = bank_htcc.getFloat("nphe",htccrow);
						//htcc_theta = bank_htcc.getFloat("theta",htccrow);		   
						//cal_energy_map = cal_cut_row(event,k);		

						if(!ec_flag){// && cal_energy_map.size() == 3){

							pcal_energy = 0;
							ecin_energy = 0;
							ecout_energy = 0;
							tot_cal_energy = 0;
							for(int ical : calMap.get(k)){
								sector = bank_cal.getByte("sector",ical);
								x_cal = bank_cal.getFloat("x",ical);
								y_cal = bank_cal.getFloat("y",ical);
								lu = bank_cal.getFloat("lu",ical);
								lv = bank_cal.getFloat("lv",ical);
								lw = bank_cal.getFloat("lw",ical);
								switch(bank_cal.getInt("layer",ical)){
									case 1:
									pcal_energy = bank_cal.getFloat("energy",ical);
									break;
									case 4:
									ecin_energy = bank_cal.getFloat("energy",ical);
									break;
									case 7:
									ecout_energy = bank_cal.getFloat("energy",ical);
									break;
									default: break;
								}
							}	


							Cal_y_vs_x_precut.fill(x_cal,y_cal);
							/*lu = cal_energy_map.get(1).get(1);
							lv = cal_energy_map.get(1).get(2);
							lw = cal_energy_map.get(1).get(3);
							pcal_energy = cal_energy_map.get(1).get(0);
							ecin_energy = cal_energy_map.get(4).get(0);
							ecout_energy = cal_energy_map.get(7).get(0);
							*/
							//System.out.println(ecin_energy + " " + ecout_energy + "      ");
							Cal_lu.fill(lu);
							Cal_lv.fill(lv);
							Cal_lw.fill(lw);				       
							if(pcal_energy != 0 && ecin_energy != 0 && ecout_energy != 0 && ecin_energy > 0.06){// && lu > 60 && lu < 350 && lv < 370 && lw < 360){
								ec_flag = true;
								tot_cal_energy = pcal_energy + ecin_energy + ecout_energy;
							}	
						}
					} 
				}
				/*	
				//dc_row = dc_cut_row(event, k);
				if(!dc_flag){
				//for(int ical : calMap.get(k)){
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

				//htcc_row = htcc_cut_row(event,k);
				if(!htcc_flag){
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
				*/
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
		/**
		*Fill Histos
		**/
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

				Vector3 n_vec_3 = new Vector3(0,0,0);
				n_vec_3.copy(p_vec_3);
				n_vec_3.negative();
				q_vec_3 = q_vec.vect();
				if(W > 1.5){
					p_hw_mom.fill(p_mom);
					p_hw_theta.fill(p_theta);
					p_hw_phi.fill(p_phi);
					xb.fill(Q2/(2*p_mass*(elec_e-elec_eprime)));
					xbc.fill(Q2/(2*((d_mass - p_vec_prime.e())*(elec_e-elec_eprime)-q_vec_3.dot(n_vec_3))));
					xbcvsx.fill(Q2/(2*p_mass*(elec_e-elec_eprime)),Q2/((2*(d_mass - p_vec_prime.e())*(elec_e-elec_eprime)-q_vec_3.dot(n_vec_3))));


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


		if(/*ec_flag &&*/p_theta > 40 && found_electron && found_proton && nphe > 2){
			ptheta_epcoin_nphe_count++;
			double eres = e_vec_prime.e()+p_vec_prime.e()-p_mass-en;

			nphevsp.fill(e_mom,nphe);
			nphevstheta.fill(e_theta,nphe);
			//Cal_lu.fill(lu);
			//Cal_lv.fill(lv);
			//Cal_lw.fill(lw);
			Cal_y_vs_x.fill(x_cal,y_cal);

			Cal_energy_p_vs_p.fill(e_mom,tot_cal_energy/e_mom);
			ECInvsOut.fill(ecout_energy,ecin_energy);
			ECvsPC.fill(pcal_energy,ecin_energy+ecout_energy);

			/*if(dc_flag && htcc_flag){ 
			DCXvsHTCCX.fill(x_dc,x_dc-x_htcc);
			DCYvsHTCCY.fill(y_dc,y_dc-y_htcc);
			DCPhivsHTCCPhi.fill(phi_dc,phi_dc-phi_htcc);
			}*/


			if(p_theta > 95){
				PhiResPEP95.fill(Math.abs(e_phi-p_phi));
				EResPEP95.fill(-eres);
				VResPEP95.fill(e_vz-p_vz);
				ZPvsZEP95.fill(e_vz,p_vz);
			}
			if(Math.abs(e_phi-p_phi) < 185 && Math.abs(e_phi-p_phi) > 175){
				phi_count++;
				EResPE_el.fill(-eres);
			}
			PhiResPE_el.fill(Math.abs(e_phi - p_phi));
			if(Math.abs(e_phi-p_phi) < 185 && Math.abs(e_phi-p_phi) > 175 && Math.abs(eres) <= 0.5){
				phi_eres_count++;
				EnPvsEnE_el.fill(e_vec_prime.e(),p_vec_prime.e());

				momentum_el.fill(e_mom);
				W_hist_el.fill(W);
				W_vs_Q2_el.fill(Q2,W);
				Phi_vs_W_el.fill(W,e_phi);

				E_vs_Theta_el.fill(e_theta,e_vec_prime.e());
				ProtonPvsTheta_el.fill(p_theta,p_mom);
				if(true){//if(W<=1){
					w_count++;
					Q2_hist_el.fill(Q2);
					//System.out.println("elastic cut!");
					for(int key : binmap.keySet()){
						if(Q2 < binmap.get(key)){binmapcount.put(key-1,binmapcount.get(key-1)+1); break;}
					}
					num_elastic_events++;
					PhiPvsPhiE_el.fill(e_phi,p_phi);

					ZPvsZE_el.fill(e_vz,p_vz);
					VResPE_el.fill(e_vz-p_vz);
					z_vs_Theta_el.fill(e_theta,e_vz);
					Phi_vs_Theta_el.fill(e_theta,e_phi);
					ThetaPvsThetaE_el.fill(e_theta,p_theta);
					ECalc_el.fill(p_mass*(1/Math.tan(p_theta*Math.PI/180)/Math.tan(e_theta*Math.PI/360) - 1));
				}
			} 
		}

		if(nphe > 2 && p_theta > 90 && p_vec_3.dot(q_vec_3) < 0 && found_proton && found_electron){
			double pq = p_vec_3.dot(q_vec_3)/(p_vec_3.mag()*q_vec_3.mag());
			LorentzVector qpmanip = new LorentzVector();
			qpmanip.copy(q_vec);
			qpmanip.sub(p_vec_prime); 
			double Wcorr = Math.sqrt((d_mass - p_vec_prime.e() + (elec_e-elec_eprime))*(d_mass - p_vec_prime.e() + (elec_e-elec_eprime))-(qpmanip.vect()).dot(qpmanip.vect()));
			double Wcorr2 = Math.sqrt((1.86 - p_vec_prime.e() + (elec_e-elec_eprime))*(d_mass - p_vec_prime.e() + (elec_e-elec_eprime))-(qpmanip.vect()).dot(qpmanip.vect()));

			if(pq < 0){
				if(pq > -0.3){
					if(p_mom > 0 && p_mom < 0.15){
						W_pbins[0].fill(W);
						Q_pbins[0].fill(Q2);
						Pprot_pbins[0].fill(p_mom);
						PprotvsThetaprot_pbins[0].fill(p_theta,p_mom);
						Wcorr_pbins[0].fill(Wcorr);
						Wcorr2_pbins[0].fill(Wcorr2);
					}else if(p_mom > 0.15 && p_mom < 0.25){
						W_pbins[1].fill(W);
						Q_pbins[1].fill(Q2);
						Pprot_pbins[1].fill(p_mom);
						PprotvsThetaprot_pbins[1].fill(p_theta,p_mom);
						Wcorr_pbins[1].fill(Wcorr);
						Wcorr2_pbins[1].fill(Wcorr2);
					}else if(p_mom > 0.25 && p_mom < 0.5){
						W_pbins[2].fill(W);
						Q_pbins[2].fill(Q2);
						Pprot_pbins[2].fill(p_mom);
						PprotvsThetaprot_pbins[2].fill(p_theta,p_mom);
						Wcorr_pbins[2].fill(Wcorr);
						Wcorr2_pbins[2].fill(Wcorr2);
					}
				}else if(pq < -0.3){
					if(p_mom > 0 && p_mom < 0.15){
						W_pbins[3].fill(W);
						Q_pbins[3].fill(Q2);
						Pprot_pbins[3].fill(p_mom);
						PprotvsThetaprot_pbins[3].fill(p_theta,p_mom);
						Wcorr_pbins[3].fill(Wcorr);
						Wcorr2_pbins[3].fill(Wcorr2);
					}else if(p_mom > 0.15 && p_mom < 0.25){
						W_pbins[4].fill(W);
						Q_pbins[4].fill(Q2);
						Pprot_pbins[4].fill(p_mom);
						PprotvsThetaprot_pbins[4].fill(p_theta,p_mom);
						Wcorr_pbins[4].fill(Wcorr);
						Wcorr2_pbins[4].fill(Wcorr2);
					}else if(p_mom > 0.25 && p_mom < 0.5){
						W_pbins[5].fill(W);
						Q_pbins[5].fill(Q2);
						Pprot_pbins[5].fill(p_mom);
						PprotvsThetaprot_pbins[5].fill(p_theta,p_mom);
						Wcorr_pbins[5].fill(Wcorr);
						Wcorr2_pbins[5].fill(Wcorr2);
					} 	     
				}
			}
		}
		found_electron = false; 
		found_proton = false;
	}
	reader.close();
}// End File Loop

boolean dc_cut(float X, float Y, int S)
{ 
	boolean result= false;
	if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
	&& (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
  
	return result;
}

/*int get_htcc_row(DataEvent event, int row){
       DataBank bank_htcc = event.getBank("REC::Cherenkov");
       int row_index = 0;
       for(int j = 0; j < bank_htcc.rows(); j++){
       	       row_index = bank_htcc.getInt("pindex",j);
	       if(row_index == row){
	       		    return j;	       
	       }
       }
       return -1;
}*/

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
Map<Integer,List<Integer>> loadMapByIndex( 
       Bank fromBank,
       String idxVarName) {
       Map<Integer,List<Integer>> map=new HashMap<Integer,List<Integer>>();
       if (fromBank!=null) {
           for (int iFrom=0; iFrom<fromBank.getRows(); iFrom++) {
               final int iTo = fromBank.getInt(idxVarName,iFrom);
               if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
               map.get(iTo).add(iFrom);
       }
    }
    return map;
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

dir.mkdir("/ep_coincidence");
dir.cd("/ep_coincidence");
dir.addDataSet(ProtonPvsTheta);
dir.addDataSet(ProtonPvsThetaPQ);
dir.addDataSet(ProtonP);
dir.addDataSet(ThetaPvsThetaE);
dir.addDataSet(ECalc);

dir.cd();
dir.mkdir("/ep_coincidence_ptheta90");
dir.cd("/ep_coincidence_ptheta90");
dir.addDataSet(ProtonP90);
dir.addDataSet(EnPvsEnE);
dir.addDataSet(ZPvsZE);
dir.addDataSet(EResPE);
dir.addDataSet(VResPE);
dir.addDataSet(momentum);
dir.addDataSet(W_vs_Q2);
dir.addDataSet(Q2_hist);
dir.addDataSet(W_hist);
dir.addDataSet(E_vs_Theta);
dir.addDataSet(z_vs_Theta);
dir.addDataSet(Phi_vs_Theta);
dir.addDataSet(Phi_vs_W);

dir.cd();
dir.mkdir("/ep_coincidence_ptheta90_W1.5");
dir.cd("/ep_coincidence_ptheta90_W1.5");
dir.addDataSet(p_hw_mom);
dir.addDataSet(p_hw_theta);
dir.addDataSet(p_hw_phi);
dir.addDataSet(xb);
dir.addDataSet(xbc);
dir.addDataSet(xbcvsx);
TGCanvas c = new TGCanvas("c","c",800,800);
c.draw(xbcvsx);
F1D func = new F1D("func","x",0,1);
c.draw(func,"same");
dir.cd();
dir.mkdir("/ep_coincidence_ptheta95_W1.5");
dir.cd("/ep_coincidence_ptheta95_W1.5");
dir.addDataSet(p_hw_momP95);
dir.addDataSet(p_hw_thetaP95);
dir.addDataSet(p_hw_phiP95);

dir.cd();
dir.mkdir("/ec_epcoincidence");
dir.cd("/ec_epcoincidence");
dir.addDataSet(PhiPvsPhiE);
dir.addDataSet(PhiResPE);

dir.cd();
dir.mkdir("/ec_epcoincidence_ptheta40");
dir.cd("/ec_epcoincidence_ptheta40");
dir.addDataSet(nphevsp);
dir.addDataSet(nphevstheta);
dir.addDataSet(Cal_energy_p_vs_p);
dir.addDataSet(ECInvsOut);
dir.addDataSet(ECvsPC);
dir.addDataSet(Cal_y_vs_x_precut);
dir.addDataSet(Cal_y_vs_x);
dir.addDataSet(Cal_lu);
dir.addDataSet(Cal_lv);
dir.addDataSet(Cal_lw);

dir.cd();
dir.mkdir("/ec_epcoincidence_ptheta95");
dir.cd("/ec_epcoincidence_ptheta95");
dir.addDataSet(EResPEP95);
dir.addDataSet(VResPEP95);
dir.addDataSet(ZPvsZEP95);
dir.addDataSet(PhiResPEP95);

dir.cd();
dir.mkdir("/ec_epcoincidence_ptheta40_phiresd5_eres<0.5");
dir.cd("/ec_epcoincidence_ptheta40_phiresd5_eres<0.5");
dir.addDataSet(W_vs_Q2_el);
dir.addDataSet(momentum_el);
dir.addDataSet(Q2_hist_el);
dir.addDataSet(W_hist_el);
dir.addDataSet(E_vs_Theta_el);
dir.addDataSet(Phi_vs_W_el);
dir.addDataSet(EnPvsEnE_el);
dir.addDataSet(PhiResPE_el);
dir.addDataSet(EResPE_el);

dir.cd();
dir.mkdir("/ec_epcoincidence_ptheta40_phiresd5_eres<0.5_W<1");
dir.cd("/ec_epcoincidence_ptheta40_phiresd5_eres<0.5_W<1");
dir.addDataSet(PhiPvsPhiE_el);
dir.addDataSet(ZPvsZE_el);
dir.addDataSet(VResPE_el);
dir.addDataSet(z_vs_Theta_el);
dir.addDataSet(Phi_vs_Theta_el);
dir.addDataSet(ThetaPvsThetaE_el);
dir.addDataSet(ECalc_el);

dir.cd();
dir.mkdir("/ProtonBins");
dir.cd("/ProtonBins");
for(int i = 0; i < 6; i++){
	dir.addDataSet(W_pbins[i]);
	dir.addDataSet(Q_pbins[i]);
	dir.addDataSet(Pprot_pbins[i]);
   	dir.addDataSet(PprotvsThetaprot_pbins[i]);
        dir.addDataSet(Wcorr_pbins[i]);
        dir.addDataSet(Wcorr2_pbins[i]);
}

dir.writeFile("testdir.hipo");
System.out.println(num_elastic_events + "/" + numevents + " Elastic events");
for(int key : binmap.keySet()){
	  float p1 = binmap.get(key);
	  float p2 = (binmap.containsKey(key+1))?binmap.get(key+1):binmap.get(key)*Math.pow(10,1/13)
	  int p3 = binmapcount.get(key);
	  System.out.println(key + " " + p1 + "-" + p2 + " " + p3);
}
System.out.println(counter + " " + banks_counter + " " + ptheta_epcoin_nphe_count + " " + phi_count + " " + phi_eres_count + " " + w_count);
System.out.println("Done!");
