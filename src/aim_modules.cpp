///////////////////////////////////////////////////////////////////////////////
//FILE: 'aim_modules.cpp'
//
//Contains all modules of class 'Aim'
//						aerodynamics()	aim[10-49]
//						propulsion()	aim[50-74]
//						seeker()		aim[75-99]
//						guidance()		aim[100-124]
//						control()		aim[125-149]
//						forces()		aim[150-159]
//						intercept()		aim[160-174]
//
// generally used variables are assigned to aim[0-9] 
//
//070412 Created by Peter H Zipfel
//130725 Building AIM5, PZi
///////////////////////////////////////////////////////////////////////////////

#include "class_hierarchy.hpp"
#include <vector>
#include <complex>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//Definition of 'aerodynamics' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[10-49]
//
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::def_aerodynamics()
{
	//Definition of module-variables
	aim[11].init("area",0,"Reference area of aim - deg","aerodynamics","data","");
	aim[14].init("alpmax",0,"Maximum angle of attack - deg","aerodynamics","data","");
	aim[20].init("alppx",0,"Total anlge of attack - deg","aerodynamics","diag","scrn,plot");
	aim[21].init("phipx",0,"Aerodynamic roll angle - deg","aerodynamics","diag","scrn,plot");
	aim[22].init("cnpaim",0,"Normal force coeff in man.plane - ND","aerodynamics","diag","");
	aim[23].init("claim",0,"Lift force coeff in velocity coor - ND","aerodynamics","out","");
	aim[24].init("cdaim",0,"Drag coeff in velocity coor - ND","aerodynamics","out","");
	aim[25].init("caaim",0,"Axial force coeff in body coor - ND","aerodynamics","out","");
	aim[26].init("cyaim",0,"Side force coeff in body coor - ND","aerodynamics","out","");
	aim[27].init("cnaim",0,"Normal force coeff in body coor - ND","aerodynamics","out","");
	aim[28].init("cnalp",0,"Normal force derivative - 1/rad","aerodynamics","out","");
	aim[29].init("cybet",0,"Side force derivative - 1/rad","aerodynamics","out","");
	aim[30].init("gmax",0,"Max g permissible, given 'alpmax' - g's","aerodynamics","out","scrn,plot");
}
///////////////////////////////////////////////////////////////////////////////
//'aerodynamics' module 
//Member function of class 'Aim'
// (1) Lift and drag coefficients from table look-up 
// (2) Converting them to body axes
// (3) Aerodynamic derivatives for autopilot
// (4) Max g's permissible
// 
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::aerodynamics()
{
	//local variables
	double phip(0);
	//local module-variables
	double alppx(0);
	double phipx(0);
	double cnpaim(0);
	double claim(0);
	double cdaim(0);
	double caaim(0);
	double cyaim(0);
	double cnaim(0);
	double cnalp(0);
	double cybet(0);
	double gmax(0);
	double cnp_max(0);
	double cdaim_max(0);
	//localizing module-variables
	//input data
	double area=aim[11].real();
	double alphax=aim[143].real();
	//input from other modules
	double grav=flat3[11].real();
	double pdynmc=flat3[13].real();
	double mach=flat3[14].real();
	int mprop=aim[50].integer();
	double mass=aim[61].real();
	double betax=aim[144].real();
	double alpmax=aim[14].real();
	//-------------------------------------------------------------------------
	//converting to aeroballistic coordinates
	double alpha=alphax*RAD;
	double beta=betax*RAD;
	double alpp=acos(cos(alpha)*cos(beta));
	double dum1=tan(beta);
	double dum2=sin(alpha);

	if(fabs(dum2)<SMALL)
		dum2=SMALL*sign(dum2);
	phip=atan2(dum1,dum2);

	//converting to degrees for output
	alppx=alpp*DEG;
	phipx=phip*DEG;

	//table look-up of lift and drag coefficient
	claim=aerotable.look_up("cl_aim_vs_alpha_mach",alppx,mach);
	if(mprop){
		cdaim=aerotable.look_up("cd_aim_on_vs_alpha_mach",alppx,mach);
	} else{
		cdaim=aerotable.look_up("cd_aim_off_vs_alpha_mach",alppx,mach);
	}
	//coefficients in body coordinates (guarding against negative values)
	double cos_alpha=cos(alpha);
	double sin_alpha=sin(alpha);
	caaim=cdaim*cos_alpha-claim*sin_alpha;
	cnpaim=cdaim*sin_alpha+claim*cos_alpha;
	cnaim=fabs(cnpaim)*cos(phip);
	cyaim=-fabs(cnpaim)*sin(phip);

	//calculating max g permissable corresponding to 'alpmax'
	double claim_max=aerotable.look_up("cl_aim_vs_alpha_mach",alpmax,mach);
	if(mprop){
		cdaim_max=aerotable.look_up("cd_aim_on_vs_alpha_mach",alpmax,mach);
	} else{
		cdaim_max=aerotable.look_up("cd_aim_off_vs_alpha_mach",alpmax,mach);
	}
	cnp_max=cdaim_max*sin(alpmax*RAD)+claim_max*cos(alpmax*RAD);

	double falphax=fabs(alphax);
	double fbetax=fabs(betax);
	if(falphax<10)
		cnalp=(0.123+0.013*falphax)*DEG;
	else
		cnalp=0.06*pow(falphax,0.625)*DEG;
	if(fbetax<10)
		cybet=-(0.123+0.013*fbetax)*DEG;
	else
		cybet=-0.06*pow(fbetax,0.625)*DEG;

	double normal_force=cnp_max*pdynmc*area;
	double weight=mass*grav;
	gmax=normal_force/weight;
	//-------------------------------------------------------------------------
	//loading module-variables
	//diagnostics
	aim[20].gets(alppx);
	aim[21].gets(phipx);
	aim[22].gets(cnpaim);
	aim[23].gets(claim);
	aim[24].gets(cdaim);
	//output to other modules
	aim[25].gets(caaim);
	aim[26].gets(cyaim);
	aim[27].gets(cnaim);
	aim[28].gets(cnalp);
	aim[29].gets(cybet);
	aim[30].gets(gmax);
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'propulsion' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[50-74]
//
// Numerical values are from AIM5
//
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::def_propulsion()
{
	aim[50].init("mprop","int",0,"Flag for propulsion modes - ND","propulsion","diag","");
	aim[51].init("pres_sl",101325,"Atmospheric pressure at SL - Pa","propulsion","data","");
	aim[52].init("aexit",0,"Nozzle exit area - m^2","propulsion","data","");
	aim[60].init("thrust",0,"Thrust at altitude - N","propulsion","out","scrn,plot");
	aim[61].init("mass",0,"Mass of missile - kg","propulsion","out","scrn,plot");
}
///////////////////////////////////////////////////////////////////////////////
//'propulsion' module
//Member function of class 'Aim'
//
// mrpop=0 motor off
//		 1 motor on
// 
// Based on AIM5
//	initial mass=63.8 kg
//
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::propulsion()
{
	//local variables

	//local module-variables
	double thrust_sl(0);
	double thrust(0);

	//localizing module-variables
	//input data
	int mprop=aim[50].integer();
	double pres_sl=aim[51].real();
	double aexit=aim[52].real();
	//input from other modules
	double time=flat3[0].real();
	double press=flat3[16].real();
	//restoring values
	double mass=aim[61].real();
	//-------------------------------------------------------------------------
	//TPo - will this work for multiple pulses?
	if(mprop==1){
		thrust_sl=proptable.look_up("thrust_vs_time",time);
		thrust=thrust_sl+(pres_sl-press)*aexit;
		mass=proptable.look_up("mass_vs_time",time);
	}
	if(time>0.0&&thrust_sl==0){
		mprop=0;
		thrust=0;
	}
	//-------------------------------------------------------------------------
	//loading module-variables
	//diagnostics
	//output to other modules
	aim[50].gets(mprop);
	aim[60].gets(thrust);
	//saving values and output
	aim[61].gets(mass);
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'seeker' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[75-99]
//
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::def_seeker()
{
	aim[1].init("acft_com_slot","int",0,"This aim slot in combus - ND","combus","out","");
	aim[2].init("VTEL",0,0,0,"Velocity of aircraft - m/s","combus","out","");
	aim[3].init("psivlx_acft",0,"Heading of aircraft - deg","combus","out","");
	aim[4].init("thtvlx_acft",0,"Flight path angle of aircraft - deg","combus","out","");
	aim[5].init("tgt_num","int",1,"Target tail # attacked by 'this' missile","combus","data","");
	aim[75].init("mseek","int",0,"Seeker: off=0, On=1 - ND","seeker","data","");
	aim[76].init("radarValid","int", 0, "Valid radar solution available", "seeker", "data", "");
	aim[78].init("last_time", 0, "Previous time step", "seeker", "data", "");
	aim[79].init("VTAEL", 0, 0, 0, "Aim aircraft wrt missile velocity - m/s", "seeker", "out", ""); 
	aim[80].init("dta",0,"Distance between aim and aircraft - m","seeker","out","scrn,plot");
	aim[81].init("dvta",0,"Closing speed of aim wrt aircraft - m/s","seeker","out","");
	aim[82].init("tgo_aim",0,"T-GO for aim to reach aircraft - s","seeker","diag","");
	aim[83].init("los_azx",0,"Az.of LOS (airc-aim) wrt aim axes - deg","seeker","diag","");
	aim[84].init("los_elx",0,"El.of LOS (airc-aim) wrt aim axes - deg","seeker","diag","");
	aim[85].init("sigdy",0,"Pitch LOS rates in aim coor - rad/s","seeker","diag","");
	aim[86].init("sigdz",0,"Yaw LOS rates in aim coor - rad/s","seeker","diag","");
	aim[87].init("UTAA",0,0,0,"Unit LOS vector in missile coord - ND","seeker","out","");
	aim[88].init("WOEA",0,0,0,"Inertial LOS rate in missile coord - rad/s","seeker","out","");
	aim[89].init("STAL",0,0,0,"Aim aircraft wrt missile position - m","seeker","out","");
	aim[90].init("radar_N", "int", 256, "Samples per radar sweep - ND", "seeker", "data", ""); // Example Definition
	aim[91].init("radar_Tm", 50e-6, "Radar sweep duration - s", "seeker", "data", "");     // Example Definition
	aim[92].init("radar_DeltaF", 75e6, "Radar frequency excursion - Hz", "seeker", "data", "");
	aim[95].init("max_power", 0, "Max RDM power", "seeker", "out", "scrn,plot");

}
///////////////////////////////////////////////////////////////////////////////
//'seeker' module
//Kinematic seeker with unlimited field-of-view
//Locked on aircraft at start of run
//Member function of class 'Aim'
//Assuming offset location in 'combus' 'Packet' at: STEL -> 2, VTEL -> 3
//
//mseek=0 seeker off
//mseek=1 seeker on
//		
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::seeker(Packet *combus,int num_vehicles,double sim_time,double int_step)
{
	//local variables
	Matrix STEL(3,1);
	Variable *data_t=NULL;
	Matrix STAL(3,1);
	Matrix VTAEL(3, 1);
	bool processingFlag = 0;
	TemporalIntegrator temporal_integrator;

	//local module-variables
	double dta(0);
	double dvta(0);
	double tgo_aim(0);
	double los_azx(0);
	double los_elx(0);
	double sigdy(0);
	double sigdz(0);
	Matrix UTAA(3,1);
	Matrix WOEA(3,1);
	int acft_com_slot(0);
	Matrix VTEL(3,1);
	double psivlx_acft(0);
	double thtvlx_acft(0);

	// Radar Signal processing variables
	double tx_gain;
	double rx_gain;
	double lambda = 0.03; // meters
	double rcs;
	double Ptx = 100.0; //watt
	double amplitude;
	double t_now = sim_time;
	double last_time = 0.0;
	double t_prev = sim_time - int_step;
	if (t_prev > 0.0) {
		last_time = aim[78].real();
		processingFlag = 1;
	}
	else {
		last_time = t_now;
	}
	bool debug_range_dump = (sim_time > 0.020 && sim_time < 0.021); // dump only at a specific time
	static double last_dump_time = -1.0;

	int		N		= aim[90].integer();
	double	T_m		= aim[91].real(); // sweep time
	double	deltaF	= aim[92].real(); // frequency excursion hz
	double	m		= deltaF / T_m;
	
	if (this->range_fft_size != N || this->range_fft_cfg == nullptr) {
		initialize_fft_resources(); // Initialize/Reinitialize if N changed or not set
	}
	//localizing module-variables
	//input data
	int tgt_num=aim[5].integer();
	int mseek=aim[75].integer();
	//input from other modules
	Matrix TAL=flat3[22].mat();
	Matrix SAEL=flat3[26].vec();
	Matrix VAEL=flat3[27].vec();
	//-------------------------------------------------------------------------
	//downloading from 'combus' target aircraft variables 
	//building aim id = t(j+1)
	int i(0);
	char number[4];
	sprintf(number,"%i",tgt_num);
	string aircraft_id="a"+string(number);
	//finding slot 'i' of aircraft in 'combus' (same as in vehicle_list)
	for(i=0;i<num_vehicles;i++)
	{
		string id=combus[i].get_id();
		if (id==aircraft_id)
		{						
			//downloading data from aircraft packet
			data_t=combus[i].get_data();
			STEL=data_t[2].vec();
			VTEL=data_t[3].vec();
			psivlx_acft=data_t[4].real();
			thtvlx_acft=data_t[5].real();
			//save aircraft com slot
			acft_com_slot=i;
		}
	}
	if(mseek){
		//LOS kinematics of target aircraft  wrt missile
		// Get previous information
		Matrix last_VTAEL = aim[79].vec();
		Matrix last_STAL = aim[89].vec();
		double last_time = aim[78].real();
		//aircraft (T) wrt missile (A) position
		STAL=STEL-SAEL;
		//unit LOS vector
		dta=STAL.absolute();
		double dum=1/dta;
		Matrix UTAL=STAL*dum;
		UTAA=TAL*UTAL;
		//LOS angles wrt missile body axes
		Matrix POLAR=UTAA.pol_from_cart();
		los_azx=DEG*POLAR[1];
		los_elx=DEG*POLAR[2];
		double los_az_rad = POLAR[1];
		double los_el_rad = POLAR[2];

		//differential velocity of target aircraft wrt aim missile as observed from Earth in local level coor
		VTAEL=VTEL-VAEL;

		// Here we have full info on the target wrt the missile
		// Insert the radar stuff here
		// dta - Unit LOS vector
		// VTAEL - differental velocity of the target
		// STAL - 3x1 position vector of the target wrt missile
		if (processingFlag) {
			// figure out the time between steps
			double delta_t = sim_time - int_step;
			// define the number of chirps that happened between the previous time and the current time
			int chirps_per_step = static_cast<int>(delta_t / T_m); // e.g., 0.002 / 0.00005 = 40


			for (int i = 0; i < chirps_per_step; ++i) {
				double ant_d = lambda / 2.0;
				double k0 = 2.0 * PI / lambda;
				double frac = static_cast<double>(i) / chirps_per_step;
				double t_chirp = last_time + frac * (t_now - last_time);

				// Interpolate STAL and VTAEL
				Matrix interp_STAL = last_STAL * (1.0 - frac) + STAL * frac;
				Matrix interp_VTAEL = last_VTAEL * (1.0 - frac) + VTAEL * frac;

				// LOS unit vector in Earth frame
				double R = interp_STAL.absolute();
				Matrix interp_UTAL = interp_STAL / R;

				// Range rate
				double Rdot = (interp_UTAL ^ interp_VTAEL);

				// Transform to missile body coordinates (IMPORTANT)
				Matrix UTAA_interp = TAL * interp_UTAL;

				// Get LOS angles for antenna/RCS lookup
				Matrix polar_interp = UTAA_interp.pol_from_cart();
				double los_az = DEG * polar_interp[1];
				double los_el = DEG * polar_interp[2];
				double d = lambda / 2.0;
				double PI = acos(-1);
				double delta_phi_az = 2 * PI * d * sin(polar_interp[1]) / lambda; // azimuth in radians
				double delta_phi_el = 2 * PI * d * sin(polar_interp[2]) / lambda; // elevation in radians
		

				// Calculate phase terms based on geometry and true angles
				double phase_az_component = (k0 * ant_d / 2.0) * sin(polar_interp[1]) * cos(polar_interp[2]);
				double phase_el_component = (k0 * ant_d / 2.0) * sin(polar_interp[2]);
				
				
				tx_gain = antenna_tx_gain(los_az, los_el);
				rx_gain = antenna_rx_gain(los_az, los_el);
				rcs = rcs_lookup(los_az, los_el);
				amplitude = radar_range_eq(Ptx, tx_gain, rx_gain, lambda, R, rcs);
				double amplitude_dbm = 10.0 * log10(amplitude * 1000.0);
				// Generate 1 sweep of 1024 samples
				ComplexVec sweep(N);
				double fb = (2 * m * R) / C;
				double Fs = N / T_m; // Sampling frequency of the ADC
				double fb_alias = fmod(fb + Fs / 2.0, Fs) - Fs / 2.0; //aliased beat freq
				double fd = (-2 * Rdot) / lambda;
				double A = sqrt(amplitude); // converting Prx to amplitude

				for (int j = 0; j < N; ++j) {
					double t = j * (T_m / N);
					double phase = 2 * PI * (fb * t - fd * t_chirp); //Non-aliased phase
					//double phase = 2 * pi * (fb_alias * t - fd * t_chirp); // aliased beat
					//sweep.push_back(std::polar(A, phase));
					sweep[j] = std::polar(A, phase); // A * exp(j*phase)
					
				}

				ComplexVec ch1 = apply_phase_offset(sweep, -phase_az_component - phase_el_component); // TR
				ComplexVec ch2 = apply_phase_offset(sweep, +phase_az_component - phase_el_component); // TL
				ComplexVec ch3 = apply_phase_offset(sweep, -phase_az_component + phase_el_component); // BR
				ComplexVec ch4 = apply_phase_offset(sweep, +phase_az_component + phase_el_component); // BL

				adc_model(ch1, 12, 1.0);
				adc_model(ch2, 12, 1.0);
				adc_model(ch3, 12, 1.0);
				adc_model(ch4, 12, 1.0);

				ComplexVec rp1 = compute_range_profile(ch1, debug_range_dump, N, T_m, m);
				ComplexVec rp2 = compute_range_profile(ch2, debug_range_dump, N, T_m, m);
				ComplexVec rp3 = compute_range_profile(ch3, debug_range_dump, N, T_m, m);
				ComplexVec rp4 = compute_range_profile(ch4, debug_range_dump, N, T_m, m);

				// Add to DSP buffer
				if (!rp1.empty()) {
					double range_resolution = C / (2 * deltaF);  // 0.3 meters // meters per bin
					double expected_bin = R / range_resolution;
					dsp_buffer4ch.add_sweep_4ch(rp1, rp2, rp3, rp4);
				}
			}
			if (dsp_buffer4ch.is_ready()) {
				DoubleMat rdm1 = dsp_buffer4ch.compute_rdm(dsp_buffer4ch.ch1); // CFAR on CH1
				DoubleMat rdm2 = dsp_buffer4ch.compute_rdm(dsp_buffer4ch.ch2);
				DoubleMat rdm3 = dsp_buffer4ch.compute_rdm(dsp_buffer4ch.ch3);
				DoubleMat rdm4 = dsp_buffer4ch.compute_rdm(dsp_buffer4ch.ch4);

				double max_power = 0;
				for (auto& row : rdm1)
					for (double val : row)
						if (val > max_power) max_power = val;
				aim[95].gets(max_power); // For diagnostics or plotting

				// CFAR Detection
				IntMat cfar_hits = log_cfar(rdm1, 2, 12, 20.0);
				// Accumulator
				temporal_integrator.add_detection_map(cfar_hits);
				IntMat integrated_detections = temporal_integrator.integrate();

				// AOA Calculations
				ComplexMat rdm_complex_1 = dsp_buffer4ch.ch1;
				ComplexMat rdm_complex_2 = dsp_buffer4ch.ch2;
				ComplexMat rdm_complex_3 = dsp_buffer4ch.ch3;
				ComplexMat rdm_complex_4 = dsp_buffer4ch.ch4;
				auto aoa_estimates = estimate_monopulse_aoa(
					rdm_complex_1,
					rdm_complex_2,
					rdm_complex_3,
					rdm_complex_4,
					integrated_detections,
					rdm1,
					lambda,
					-60.0,
					los_az_rad,
					los_el_rad);

				AoA_Result bestDetection = get_best_detection(aoa_estimates, rdm1);
				if (bestDetection.valid) {
					const AoA_Estimate& det = bestDetection.estimate;

					// Convert angle estimates from radians to degrees
					double est_az_deg = DEG * det.azimuth_est;
					double est_el_deg = DEG * det.elevation_est;

					// Estimate physical range and velocity
					double range_resolution = C / (2 * deltaF);
					double CPI = 32 * T_m; // CPI = max_sweeps * T_m
					double doppler_resolution = 1/CPI; //Hz
					double R_est = det.range_bin * range_resolution;
					// Velocity wrapping
					int N_doppler = 32; // Should use max_sweeps
					double fd_est = 0.0;
					if (det.doppler_bin < N_doppler / 2) { // Positive frequencies
						fd_est = det.doppler_bin * doppler_resolution;
					}
					else { // Negative frequencies (wrapped)
						fd_est = (det.doppler_bin - N_doppler) * doppler_resolution;
					}
					double V_est = fd_est * lambda / 2.0; // Positive V = Towards radar
					int radarSolution = 1;
					// Save to aim variables for downstream use (replace with correct index if needed)
					aim[76].gets(radarSolution);		  // Valid Radar Solution
					
					// Estimate LOS Vector in body frame from AOA
					Matrix UTAA_est(3, 1);
					// Need to flip the elevation to get the proper utaa
					UTAA_est.cart_from_pol(1.0, det.azimuth_est, -det.elevation_est); // <<< Potential issue (2 - Sign Convention)

					// Estimate target Position in body frame using R_est and UTAA
					Matrix STAL_est = UTAA_est*R_est;

					// Estimate relative VELOCITY in BODY frame (VTAA_est)
					// We have V_est (scalar velocity along LOS) and UTAA_est (LOS direction in body)
					// Assumption: Treat estimated velocity V_est as the velocity along the estimated LOS UTAA_est.
					// VTAA = Velocity of Target relative to AIM, in AIM body frame
					Matrix VTAA_est = UTAA_est * V_est; // <<< POTENTIAL ISSUE 3 (Only LOS component)

					// Convert estimated relative velocity FROM BODY frame TO EARTH frame (VTAEL_est)
					// Need T.M. from Body to Earth, which is TAL.trans()
					// VTAEL = Velocity of Target relative to AIM, in EARTH frame
					Matrix VTAEL_est = TAL.trans() * VTAA_est;

					// Convert UTAA from body to earth frame
					Matrix UTAL_est = TAL.trans() * UTAA_est;

					// Calculate closing velocity and tgo
					double dvta_est = UTAL_est ^ VTAEL_est;
					//dvta_est = V_est; // More direct
					double tgo_aim_est = R_est / fabs(dvta_est);

					// Calculate inertial LOS rate in body frame
					Matrix WOEA_est = TAL * UTAL_est.skew_sym() * VTAEL_est * (1.0 / R_est);

					
					aim[96].gets(R_est);      // Estimated range
					aim[97].gets(V_est);      // Estimated closing velocity
					aim[98].gets(est_az_deg); // Estimated azimuth
					aim[99].gets(est_el_deg); // Estimated elevation
					aim[81].gets(dvta_est);			// check for radar solution
					aim[87].gets_vec(UTAA_est);		// check for radar solution
					//aim[88].gets_vec(WOEA_est);		// check for radar solution

					// Optionally, you could build an estimated LOS vector using polar coords here
				}
				else {
					std::cout << "No AOA Estimates" << std::endl;
				}
				//std::cout << "Filtered AoA detections: " << aoa_estimates.size() << std::endl;
				// Optionally dump or plot
				
				if (sim_time - last_dump_time >= 0.1) {
					dump_rdm_power(rdm1, "rdm_power.dat");
					dump_rdm_binary(cfar_hits, "cfar_binary.dat");
					dump_rdm_binary(integrated_detections, "rdm_integrated.dat");
					dump_aoa_estimates(aoa_estimates, "aoa_estimates.dat");
					last_dump_time = sim_time;
				}
				
				
			}

		}


		//closing velocity
		dvta = UTAL ^ VTAEL;
		//time-to-go
		tgo_aim = dta / fabs(dvta);

		//inertial LOS rates in missile body coordinates
		WOEA = TAL * UTAL.skew_sym() * VTAEL * dum;
		//diagnostic
		sigdy = WOEA[1];
		sigdz = WOEA[2];
	}
	//-------------------------------------------------------------------------
	//loading module-variables
	//output to other modules
	int radarSolution = aim[76].integer();
	if (!radarSolution) {
		aim[81].gets(dvta);			// check for radar solution
		aim[87].gets_vec(UTAA);		// check for radar solution
		//aim[88].gets_vec(WOEA);		// check for radar solution
	}
	aim[88].gets_vec(WOEA);		// check for radar solution

	aim[1].gets(acft_com_slot); // no change
	aim[2].gets_vec(VTEL); // no change
	aim[3].gets(psivlx_acft); // no change
	aim[4].gets(thtvlx_acft); // no change
	aim[89].gets_vec(STAL);		// check for radar solution
	aim[78].gets(last_time);	// no change
	aim[79].gets_vec(VTAEL);	// no change
	//saving values and output
	aim[83].gets(los_azx);		// check for radar solution
	aim[84].gets(los_elx);		// check for radar solution
	//diagnostics
	aim[80].gets(dta);			// check for radar solution
	aim[82].gets(tgo_aim);		// check for radar solution
	aim[85].gets(sigdy);		// check for radar solution
	aim[86].gets(sigdz);		// check for radar solution
	
}
///////////////////////////////////////////////////////////////////////////////
//'radar Processing' support functions
//Kinematic seeker with unlimited field-of-view
//Locked on aircraft at start of run
//Member function of class 'Aim'
//
//		
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////
double Aim::antenna_tx_gain(double az, double el) 
{
	// Simple model: fixed gain or bilinear interpolation from table
	double gainDb = gaintable.look_up("az_vs_el_db", az, el);
	return gainDb; // dB or linear gain
}
double Aim::antenna_rx_gain(double az, double el)
{
	// Simple model: fixed gain or bilinear interpolation from table
	double gainDb = gaintable.look_up("az_vs_el_db", az, el);
	return gainDb; // dB or linear gain
}
double Aim::rcs_lookup(double az, double el) {
	// Convert to appropriate frame, lookup in RCS table, interpolate
	return 0.1; // m^2
}
double Aim::radar_range_eq(double Ptx, double Gtx, double Grx, double lambda, double R, double sigma) {
	double pi = acos(-1);
	double num = Ptx * Gtx * Grx * pow(lambda, 2) * sigma;
	double denom = pow(4 * pi, 3) * pow(R, 4);
	return num / denom; // received power in Watts
}

void Aim::synth_signal(double t_now, double m, double lambda, double Tm, double R, double Rdot, double A)
{
	double c = 3E8;
	double pi = acos(-1);
	double fb = (2 * m * R) / c;
	double fd = (-2 * Rdot) / lambda;
	double t_chirp = fmod(t_now, Tm);
	double phase = 2 * pi * (fb * t_chirp - fd * t_now);
	double real_sig = A * cos(phase);
	double imag_sig = A * sin(phase);

	// Optionally: store in aim[]
	aim[90].gets(real_sig);
	aim[91].gets(imag_sig);
	aim[92].gets(fb);
	aim[93].gets(fd);
}

AoA_Result Aim::get_best_detection(const std::vector<AoA_Estimate>& estimates, const DoubleMat& rdm_power) {
	AoA_Result result;
	result.valid = false;

	double max_power = -1.0;
	AoA_Estimate best;
	const AoA_Estimate* best_ptr = nullptr;
	for (const auto& est : estimates) {
		// Bounds checking is important!
		if (est.range_bin >= rdm_power.size() || est.doppler_bin >= rdm_power[est.range_bin].size()) {
			std::cerr << "*** Warning: Index out of bounds in get_best_detection for range "
				<< est.range_bin << ", doppler " << est.doppler_bin << " ***" << std::endl;
			continue; // Skip this invalid estimate index
		}
		double power = rdm_power[est.range_bin][est.doppler_bin];

		if (power > max_power) {
			max_power = power;
			// best = est; // Instead of copying here...
			best_ptr = &est; // ...store a pointer to the best estimate found
			result.valid = true;
		}
	}
	// --- Assign ONLY if a valid estimate was found ---
	if (result.valid) {
		// A valid estimate was found, assign it from the pointer
		result.estimate = *best_ptr; // Dereference the pointer to copy the best estimate
	}
	else {
		// No valid estimate found (estimates was empty or all below threshold)
		// result.estimate remains uninitialized, but result.valid is false,
		// so the caller should check result.valid before using result.estimate.
		// Optionally, initialize result.estimate to some default here if desired:
		// result.estimate = AoA_Estimate(); // Requires default constructor for AoA_Estimate
	}
	return result;
}
///////////////////////////////////////////////////////////////////////////////
//'adc_model'  function
//Applies an AGC and quantizes the signal.
//Moves from real, physical power to dBFS
//		
//070412 Created by Peter H Zipfel
//130727 Modified for AIM5, PZi
///////////////////////////////////////////////////////////////////////////////

void Aim::adc_model(ComplexVec& sweep, int bit_depth, double voltage_range)
{
	const size_t N = sweep.size();
	const double half_range = voltage_range * 0.5;
	double max_magnitude_sq = 0.0;

	// Step 1: AGC peak search using squared magnitude (avoids sqrt)
	for (const auto& s : sweep) {
		double mag_sq = s.real() * s.real() + s.imag() * s.imag();
		if (mag_sq > max_magnitude_sq)
			max_magnitude_sq = mag_sq;
	}

	// Compute AGC gain
	double gain = 1.0;
	if (max_magnitude_sq > 0.0) {
		double peak = std::sqrt(max_magnitude_sq);
		gain = half_range / peak;
	}

	// Step 2: Quantize
	const int max_code = (1 << (bit_depth - 1)) - 1;
	const int min_code = -max_code;
	const double scale_factor = max_code / half_range;

	for (size_t i = 0; i < N; ++i) {
		double real = sweep[i].real() * gain;
		double imag = sweep[i].imag() * gain;

		int q_real = static_cast<int>(real * scale_factor + 0.5);
		int q_imag = static_cast<int>(imag * scale_factor + 0.5);

		// Clamp
		if (q_real > max_code) q_real = max_code;
		if (q_real < min_code) q_real = min_code;
		if (q_imag > max_code) q_imag = max_code;
		if (q_imag < min_code) q_imag = min_code;

		sweep[i] = Complex(q_real / scale_factor, q_imag / scale_factor);
	}
}
///////////////////////////////////////////////////////////////////////////////
//'match_filter'  function
//Standard match filtering function
//Maximizes SNR for a known signal (up-chirp)
//		Convolution with the time-reversed complex conjugate of the transmitted
//		signal
//250412 Savas N Mavridis
///////////////////////////////////////////////////////////////////////////////

ComplexVec Aim::matched_filter(const ComplexVec& sweep, const ComplexVec& template_chirp){
	int N = sweep.size();
	std::vector<std::complex<double>> output(N, { 0, 0 });

	for (int n = 0; n < N; ++n) {
		std::complex<double> acc(0, 0);
		for (int k = 0; k <= n; ++k) {
			acc += sweep[k] * template_chirp[n - k];
		}
		output[n] = acc;
	}

	return output;
}

///////////////////////////////////////////////////////////////////////////////
//'match_filter_template'  function
//Match Filter Template
//
//250412 Savas N Mavridis
///////////////////////////////////////////////////////////////////////////////
ComplexVec Aim::generate_chirp_template(int N, double Tm, double m) {
	ComplexVec template_chirp(N); // Pre-allocate size
	const double PI = acos(-1.0);

	for (int i = 0; i < N; ++i) {
		double t = i * (Tm / N);
		double phase = PI * m * t * t; // Phase = pi * m * t^2
		template_chirp[i] = std::exp(std::complex<double>(0.0, phase)); // exp(j * phase)
	}
	return template_chirp;
}

///////////////////////////////////////////////////////////////////////////////
// 'compute_range_profile' function (Revised Style)
// 
// Performs dechirp, windowing, and FFT on a single sweep.
// Receives all parameters as arguments.
// 
///////////////////////////////////////////////////////////////////////////////
/* 
ComplexVec Aim::compute_range_profile(const ComplexVec& sweep, bool debug_dump, int N, double Tm, double m)
{
	// --- Use Passed-in Parameters ---
	// N, Tm, m are now function arguments
	const double PI = acos(-1.0);
	// Input Validation
	if (sweep.size() != N) {
		std::cerr << "*** Error: Sweep size (" << sweep.size()
			<< ") does not match N (" << N << ") in compute_range_profile ***" << std::endl;
		return ComplexVec();
	}

	// --- Windowing (Apply Hann Window) ---
	ComplexVec windowed_sweep(N);
	for (int k = 0; k < N; ++k) {
		double hann_mult = 0.5 * (1.0 - cos(2.0 * PI * k / (N - 1)));
		windowed_sweep[k] = sweep[k] * hann_mult;
	}

	// --- FFT ---
	ComplexVec range_fft = kiss_fft_wrapper(windowed_sweep); // Assumes fft is available

	// --- Optional Debug Dump ---
	if (debug_dump) {
		dump_range_profile(range_fft, "range_profile_dechirp.dat"); // Assumes dump_range_profile is available
	}

	return range_fft;
}
*/
ComplexVec Aim::compute_range_profile(const ComplexVec& sweep, bool debug_dump, int N, double Tm, double m)
{
	// ... check input size N against this->range_fft_size ...
	if (N != this->range_fft_size || !this->range_fft_cfg) {
		std::cerr << "*** Error: Range FFT resources not initialized or size mismatch in compute_range_profile ***" << std::endl;
		return ComplexVec();
	}

	const double PI = acos(-1.0);
	ComplexVec windowed_sweep(N);
	for (int k = 0; k < N; ++k) {
		double hann_mult = 0.5 * (1.0 - cos(2.0 * PI * k / (N - 1)));
		windowed_sweep[k] = sweep[k] * hann_mult;
	}

	// --- FFT ---
	// Use the new execution function with pre-allocated resources
	ComplexVec range_fft(N); // Vector to store results
	bool success = kiss_fft_exec(this->range_fft_cfg, this->range_fft_in, this->range_fft_out, N, windowed_sweep);

	if (success) {
		// Copy result from kiss_fft_cpx* out buffer to ComplexVec result
		for (int k = 0; k < N; ++k) {
			range_fft[k] = Complex(this->range_fft_out[k].r, this->range_fft_out[k].i);
		}
	}
	else {
		// Handle FFT execution error, maybe return empty vector
		return ComplexVec();
	}

	// --- Optional Debug Dump ---
	if (debug_dump) {
		dump_range_profile(range_fft, "range_profile_dechirp.dat");
	}

	return range_fft;
}


///////////////////////////////////////////////////////////////////////////////
// 'ca_cfar' function 
// 
// Simple cell averaging constant false alarm rate detector
// 
///////////////////////////////////////////////////////////////////////////////
IntMat Aim::ca_cfar(const DoubleMat& rdm, int guard_cells = 2, int training_cells = 8, double threshold_scale = 5.0){  // scaling factor over noise
	int rows = rdm.size();
	int cols = rdm[0].size();
	IntMat detections(rows, std::vector<int>(cols, 0));

	for (int i = training_cells + guard_cells; i < rows - training_cells - guard_cells; ++i) {
		for (int j = training_cells + guard_cells; j < cols - training_cells - guard_cells; ++j) {
			double noise_sum = 0;
			int count = 0;

			for (int m = i - training_cells - guard_cells; m < i - guard_cells; ++m)
				for (int n = j - training_cells - guard_cells; n < j + guard_cells + 1; ++n)
					if (m != i || n != j) {
						noise_sum += rdm[m][n];
						count++;
					}

			for (int m = i + guard_cells + 1; m <= i + training_cells + guard_cells; ++m)
				for (int n = j - training_cells - guard_cells; n < j + guard_cells + 1; ++n)
					if (m != i || n != j) {
						noise_sum += rdm[m][n];
						count++;
					}

			double noise_mean = noise_sum / count;
			double threshold = noise_mean * threshold_scale;

			if (rdm[i][j] > threshold)
				detections[i][j] = 1;
		}
	}

	return detections;
}
///////////////////////////////////////////////////////////////////////////////
// 'log_cfar' function 
// 
// Logarithmic cell averaging constant false alarm rate detector
// 
///////////////////////////////////////////////////////////////////////////////

IntMat Aim::log_cfar(const DoubleMat& rdm, int guard_cells, int training_cells, double threshold_dB)
{
	int rows = rdm.size();
	int cols = rdm[0].size();
	IntMat detections(rows, IntVec(cols, 0));

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {

			// Check bounds
			if (i < training_cells + guard_cells || i >= rows - (training_cells + guard_cells))
				continue;
			if (j < training_cells + guard_cells || j >= cols - (training_cells + guard_cells))
				continue;

			// Collect training cells (in dB)
			std::vector<double> training_db;
			for (int m = i - training_cells - guard_cells; m <= i + training_cells + guard_cells; ++m) {
				for (int n = j - training_cells - guard_cells; n <= j + training_cells + guard_cells; ++n) {
					// Skip guard cells and CUT
					if (std::abs(m - i) <= guard_cells && std::abs(n - j) <= guard_cells)
						continue;
					if (m == i && n == j)
						continue;
					double val = rdm[m][n];
					if (val > 0.0)
						training_db.push_back(10.0 * std::log10(val));
				}
			}

			if (training_db.empty()) continue;

			// Compute noise estimate in dB
			double sum_db = 0.0;
			for (double db : training_db)
				sum_db += db;
			double mean_db = sum_db / training_db.size();

			// Convert CUT to dB
			double cut_db = 10.0 * std::log10(rdm[i][j]);

			// Detection if CUT exceeds threshold
			if (cut_db > mean_db + threshold_dB)
				detections[i][j] = 1;
		}
	}

	return detections;
}

///////////////////////////////////////////////////////////////////////////////
// 'apply_phase_offset' function 
// 
// Applies directional phase offset
// 
///////////////////////////////////////////////////////////////////////////////
/*
ComplexVec Aim::apply_phase_offset(const ComplexVec& sweep, double phase_offset_rad) {
	ComplexVec shifted(sweep.size());
	for (size_t k = 0; k < sweep.size(); ++k) {
		shifted[k] = sweep[k] * std::exp(Complex(0, phase_offset_rad)); // rotate by j \delta \phi
	}
	return shifted;
}
*/

ComplexVec Aim::apply_phase_offset(const ComplexVec& input, double delta_phi) {
	ComplexVec output = input;
	Complex rotator = std::polar(1.0, delta_phi);
	for (auto& s : output)
		s *= rotator;
	return output;
}



void Aim::initialize_fft_resources() {
	// Clean up existing resources first, if any (e.g., if re-initializing)
	if (range_fft_cfg) free(range_fft_cfg);
	delete[] range_fft_in;
	delete[] range_fft_out;
	range_fft_cfg = nullptr; // Reset pointers after deleting
	range_fft_in = nullptr;
	range_fft_out = nullptr;

	// Get N from the aim array - CRITICAL: Ensure aim[90] is valid!
	// Add error handling for accessing aim[90] if needed
	if (90 >= NAIM) { // Basic bounds check assuming NAIM is max index+1
		std::cerr << "*** Error: Index 90 out of bounds for aim array in initialize_fft_resources ***" << std::endl;
		range_fft_size = 0;
		return;
	}
	range_fft_size = aim[90].integer();

	if (range_fft_size <= 0) {
		std::cerr << "*** Error: Invalid range FFT size N = " << range_fft_size << " read from aim[90] ***" << std::endl;
		// range_fft_cfg, _in, _out are already nullptr
		return; // Prevent allocation with bad size
	}

	std::cout << "Initializing Range FFT resources for size N = " << range_fft_size << std::endl;
	try {
		range_fft_in = new kiss_fft_cpx[range_fft_size];
		range_fft_out = new kiss_fft_cpx[range_fft_size];
		range_fft_cfg = kiss_fft_alloc(range_fft_size, 0, nullptr, nullptr); // 0 for forward FFT

		if (!range_fft_cfg) { // kiss_fft_alloc returns null on failure
			throw std::bad_alloc(); // Treat allocation failure consistently
		}
	}
	catch (const std::bad_alloc& e) {
		std::cerr << "*** Error: Failed to allocate KissFFT resources for Range FFT: " << e.what() << " ***" << std::endl;
		// Clean up partially allocated resources
		delete[] range_fft_in; // delete[] on nullptr is safe
		delete[] range_fft_out;
		if (range_fft_cfg) free(range_fft_cfg); // free(nullptr) might be unsafe, check if needed
		range_fft_cfg = nullptr; range_fft_in = nullptr; range_fft_out = nullptr;
		range_fft_size = 0; // Indicate failure
		// Consider re-throwing or handling error more formally
	}
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'guidance' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[100-124]
//
//	mguid = |manvr|mode|
//			  manvr = 1 decaying spiral maneuver	
//				    mode =1 pro-nav 
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::def_guidance()
{
	//Definition of module-variables
	aim[100].init("mguid","int",0,"=|manvr|mode| =11:|spiral|pronav|","guidance","data","");
	aim[101].init("gnav",0,"Proportional navigation gain - ND","guidance","data","");
	aim[102].init("tgo_manvr",0,"Evasive maneuver start - sec","guidance","data","");
	aim[103].init("amp_manvr",0,"Initial amplitude of maneuver - g's","guidance","data","");
	aim[104].init("frq_manvr",0,"Frequency of spiral maneuver - rad/s","guidance","data","");
	aim[105].init("tgo63_manvr",0,"Amplitude decay at 63% tgo - sec","guidance","data","");
	aim[106].init("annx",0,"Normal accel command, unrestricted - g's","guidance","diag","");
	aim[107].init("allx",0,"Lateral accel command, unrestricted  - g's","guidance","diag","");
	aim[108].init("an_manvr",0,"Normal spiral command  - g's","guidance","diag","");
	aim[109].init("al_manvr",0,"Lateral spiral command  - g's","guidance","diag","");
	aim[110].init("ancomx",0,"Aim normal acceleration command - g's","guidance","out","scrn,plot");
	aim[111].init("alcomx",0,"Aim lateral acceleration command - g's","guidance","out","scrn,plot");
}
///////////////////////////////////////////////////////////////////////////////
//'guidance' module 
//Member function of class 'Aim' 
//Calculating commands for autopilot
//
//Spiral decay maneuver
//  Starts at 'tgo_manvr' with amplitude 'amp_manvr' and frequency 'frq_manvr'.
//	Amplitude decays exponentially with the 63% time-to-go value 'tgo63_manvr' 
//
//070412 Created by Peter H Zipfel
//070905 Added spiral evasive maneuver, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::guidance(Packet *combus,int num_vehicles)
{
	static bool header_flag = true;
	//local variables
	double phi(0);

	//local module-variables
	double ancomx(0);
	double alcomx(0);
	double annx(0);
	double allx(0);
	double amp(0);
	double an_manvr(0);
	double al_manvr(0);

	//localizing module-variables
	//input data
	int mguid=aim[100].integer();
	double gnav=aim[101].real();
	double tgo_manvr=aim[102].real();
	double amp_manvr=aim[103].real();
	double frq_manvr=aim[104].real();
	double tgo63_manvr=aim[105].real();
	//input from other modules
	double sim_time = flat3[0].real(); // Get simulation time for logging
	double grav=flat3[11].real();
	double gmax=aim[30].real();
	double dvta=aim[81].real();
	double tgo_aim=aim[82].real();
	Matrix UTAA=aim[87].vec();
	Matrix WOEA=aim[88].vec();
	//-------------------------------------------------------------------------
	//decoding guidance flag
	int guid_manvr=mguid/10;
	int guid_mode=(mguid%10);
	
	// Log input data
	try {
		std::ofstream log_file;
		log_file.open("guidanceLog.txt", std::ios::app);

		if (log_file.is_open()) 
		{
			if (header_flag) 
			{
				log_file << "simTime\tGrav\tGmax\tdvta\ttgo_aim\t"
					<< "UTAA_X\tUTAA_Y\tUTAA_Z\t"
					<< "WOEA_X\tWOEA_Y\tWOEA_Z\n";
				header_flag = false;
			}
			log_file << std::fixed << std::setprecision(6);
			// Write current data values, separated by tabs
			log_file << sim_time << "\t"
				<< grav << "\t"
				<< gmax << "\t"
				<< dvta << "\t"
				<< tgo_aim << "\t"
				// Extract components from Matrix objects (assuming 3x1)
				<< UTAA.get_loc(0, 0) << "\t"
				<< UTAA.get_loc(1, 0) << "\t"
				<< UTAA.get_loc(2, 0) << "\t"
				<< WOEA.get_loc(0, 0) << "\t"
				<< WOEA.get_loc(1, 0) << "\t"
				<< WOEA.get_loc(2, 0) << "\n";

			log_file.close(); // Close the file after writing

		}
		else 
		{
			// Optional: Print error message (maybe only once)
			static bool log_error_shown = false;
			if (!log_error_shown) {
				std::cerr << "Warning: Could not open " << "guidanceLog.txt" << " for logging guidance inputs." << std::endl;
				log_error_shown = true;
			}

		}
	}
	catch (const std::exception& e) {
		// Optional: Catch potential file exceptions
		static bool log_exception_shown = false;
		if (!log_exception_shown) {
			std::cerr << "Exception while writing to guidance log: " << e.what() << std::endl;
			log_exception_shown = true;
		}
	}

	//pro-nav guidance
	if(guid_mode==1){
		//acceleration command in missile axes
		Matrix APNA=WOEA.skew_sym()*UTAA*gnav*fabs(dvta);
		annx=-APNA.get_loc(2,0)/grav;
		allx=APNA.get_loc(1,0)/grav;
	}
	//adding exponentially decaying spiral maneuver
	if(guid_manvr==1&&tgo_aim<tgo_manvr){
		amp=amp_manvr*(1-exp(-tgo_aim/tgo63_manvr));
		an_manvr=amp*sin(frq_manvr*tgo_aim);
		al_manvr=amp*cos(frq_manvr*tgo_aim);
		annx+=an_manvr;
		allx+=al_manvr;
	}
	//limiting acceleration commands by circular limiter
	double aax=sqrt(allx*allx+annx*annx);
	if(aax>gmax) aax=gmax;
	if((fabs(annx)<SMALL||fabs(allx)<SMALL))
		phi=0;
	else{
		phi=atan2(annx,allx);
	}
	alcomx=aax*cos(phi);
	ancomx=aax*sin(phi);
	//-------------------------------------------------------------------------
	//loading module-variables
	//output to other modules
	aim[106].gets(annx);
	aim[107].gets(allx);
	aim[108].gets(amp);
	aim[110].gets(ancomx);
	aim[111].gets(alcomx);
	//diagnostics
	aim[106].gets(annx);
	aim[107].gets(allx);
	aim[108].gets(an_manvr);
	aim[109].gets(al_manvr);
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'control' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[125-149]
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::def_control()
{
	//Definition of module-variables
	aim[127].init("ta",0,"Ratio of prop/integral gain - ND","control","data","");
	aim[128].init("tr",0,"Rate loop time constant - sec","control","data","");
	aim[129].init("gacp",0,"Root locus gain of accel loop - rad/s2","control","data","");
	aim[130].init("tip",0,"Incidence lag time constant of aim - sec","control","diag","scrn,plot");
	aim[131].init("xi",0,"Integral feedback  - rad/s","control","state","");
	aim[132].init("xid",0,"Integral feedback derivative - rad/s^2 ","control","state","");
	aim[133].init("ratep",0,"Pitch rate  - rad/s","control","state","");
	aim[134].init("ratepd",0,"Pitch rate derivative  - rad/s^2","control","state","");
	aim[135].init("alp",0,"Angle of attack - rad","control","state","");
	aim[136].init("alpd",0,"Angle of attack derivative - rad/s","control","state","");
	aim[137].init("yi",0,"Integral feedback  - rad/s","control","state","");
	aim[138].init("yid",0,"Integral feedback derivative - rad/s^2 ","control","state","");
	aim[139].init("ratey",0,"Yaw rate  - rad/s","control","state","");
	aim[140].init("rateyd",0,"Yaw rate derivative  - rad/s^2","control","state","");
	aim[141].init("bet",0,"Sideslip angle - rad","control","state","");
	aim[142].init("betd",0,"Sideslip angle - rad/s","control","state","");
	aim[143].init("alphax",0,"Angle of attack of aim - deg","control","in/out","scrn,plot");
	aim[144].init("betax",0,"Sideslip angle of aim - deg","control","in/out","scrn.plot");
}
///////////////////////////////////////////////////////////////////////////////
//Initial calculations of 'control' module 
//Member function of class 'Aim'
// 
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::init_control()
{
	//localizing module-variables
	//input data
	double alphax=aim[143].real();
	double betax=aim[144].real();
	//-------------------------------------------------------------------------
	// initializing incidence angles
	double alp=alphax*RAD;
	double bet=betax*RAD;
	//-------------------------------------------------------------------------
	//initializing variables of 'control' module
	aim[135].gets(alp);
	aim[141].gets(bet);
}///////////////////////////////////////////////////////////////////////////////
//'control' module 
//Member function of class 'Aim' 
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::control(double int_step)
{
	//local module-variables
	double tip(0);
	double alphax(0);
	double betax(0);
	//localizing module-variables
	//input data
	double ta=aim[127].real();
	double tr=aim[128].real();
	double gacp=aim[129].real();
	//input from other modules
	double grav=flat3[11].real();
	double pdynmc=flat3[13].real();
	double dvae=flat3[25].real();
	double area=aim[11].real();
	double alpmax=aim[14].real();
	double cyaim=aim[26].real();
	double cnaim=aim[27].real();
	double cnalp=aim[28].real();
	double cybet=aim[29].real();
	double thrust=aim[60].real();
	double mass=aim[61].real();
	double ancomx=aim[110].real();
	double alcomx=aim[111].real();
	//state variables
	double xi=aim[131].real();
	double xid=aim[132].real();
	double ratep=aim[133].real();
	double ratepd=aim[134].real();
	double alp=aim[135].real();
	double alpd=aim[136].real();
	double yi=aim[137].real();
	double yid=aim[138].real();
	double ratey=aim[139].real();
	double rateyd=aim[140].real();
	double bet=aim[141].real();
	double betd=aim[142].real();
	//-------------------------------------------------------------------------
	//Pitch acceleration controller
	//incidence lag time constant
	tip=dvae*mass/(pdynmc*area*fabs(cnalp)+thrust);
	//pitch specific force
	double fspz=-pdynmc*area*cnaim/mass;
	//P-I shaping
	double gr=gacp*tip*tr/dvae;
	double gi=gr/ta;
	double abez=-ancomx*grav;
	double ep=abez-fspz;
	double xid_new=gi*ep;
	xi=integrate(xid_new,xid,xi,int_step);
	xid=xid_new;
	double ratepc=-(ep*gr+xi);
	//pitch rate first order lag
	double ratepd_new=(ratepc-ratep)/tr;
	ratep=integrate(ratepd_new,ratepd,ratep,int_step);
	ratepd=ratepd_new;
	//incidence lag
	double alpd_new=(tip*ratep-alp)/tip;
	alp=integrate(alpd_new,alpd,alp,int_step);
	alpd=alpd_new;
	alphax=alp*DEG;
	//alpha limiter
	if(fabs(alphax)>alpmax) alphax=alpmax*sign(alphax);

	//Yaw acceleration controller
	//incidence lag time constant
	double tiy=dvae*mass/(pdynmc*area*fabs(cybet)+thrust);
	//yaw specific force
	double fspy=pdynmc*area*cyaim/mass;
	//P-I shaping
	gr=gacp*tiy*tr/dvae;
	gi=gr/ta;
	double abey=alcomx*grav;
	double ey=abey-fspy;
	double yid_new=gi*ey;
	yi=integrate(yid_new,yid,yi,int_step);
	yid=yid_new;
	double rateyc=(ey*gr+yi);
	//yaw rate first order lag
	double rateyd_new=(rateyc-ratey)/tr;
	ratey=integrate(rateyd_new,rateyd,ratey,int_step);
	rateyd=rateyd_new;
	//incidence lag
	double betd_new=-(tiy*ratey+bet)/tiy;
	bet=integrate(betd_new,betd,bet,int_step);
	betd=betd_new;
	betax=bet*DEG;
	//beta limiter (same as alpha limiter)
	if(fabs(betax)>alpmax) betax=alpmax*sign(betax);
	//-------------------------------------------------------------------------
	//loading module-variables
	//state variables
	aim[131].gets(xi);
	aim[132].gets(xid);
	aim[133].gets(ratep);
	aim[134].gets(ratepd);
	aim[135].gets(alp);
	aim[136].gets(alpd);
	aim[137].gets(yi);
	aim[138].gets(yid);
	aim[139].gets(ratey);
	aim[140].gets(rateyd);
	aim[141].gets(bet);
	aim[142].gets(betd);
	//output to other modules
	aim[143].gets(alphax);
	aim[144].gets(betax);
	//diagnostics
	aim[130].gets(tip);
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'force' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[150-159]
//
//Note that FSPA is entered into the 'flat3[20]' array because it is needed
// for the 'newton' module, which is a member of the 'Flat3' class
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::def_forces()
{
	//Definition of module-variables
	flat3[20].init("FSPA",0,0,0,"Specific force in vehicle coor - m/s^2","forces","out","");
	aim[150].init("aax",0,"Axial acceleration of aim - g's","forces","diag","");
	aim[151].init("alx",0,"Yaw maneuver acceleration of aim - g's","forces","diag","scrn,plot");
	aim[152].init("anx",0,"Pitch maneuver acceleration of aim - g's","forces","diag","scrn,plot");
}
///////////////////////////////////////////////////////////////////////////////
//'force' module 
//Member function of class 'Aim' 
//Calculates forces acting on the missile
//
//070412 Created by Peter H Zipfel
//080606 Added fixed aim option, which requries g-bias in vertical direction, PZi
///////////////////////////////////////////////////////////////////////////////
void Aim::forces()
{
	//local module-variables
	Matrix FSPA(3,1);
	double aax(0);
	double alx(0);
	double anx(0);

	//localizing module-variables
	//input data
	double acc_longx=aim[40].real();
	//input from other modules
	double grav=flat3[11].real();
	double pdynmc=flat3[13].real();
	double area=aim[11].real();
	double caaim=aim[25].real();
	double cyaim=aim[26].real();
	double cnaim=aim[27].real();
	double thrust=aim[60].real();
	double mass=aim[61].real();
	//-------------------------------------------------------------------------
	FSPA[0]=(thrust-caaim*pdynmc*area)/mass;
	FSPA[1]=(cyaim*pdynmc*area)/mass;
	FSPA[2]=(-cnaim*pdynmc*area)/mass;

	//diagnostics: accelerations in aim body coord
	aax=FSPA[0]/grav;
	alx=FSPA[1]/grav;
	anx=-FSPA[2]/grav;
	//-------------------------------------------------------------------------
	//loading module-variables
	//output to other modules
	flat3[20].gets_vec(FSPA);
	//diagnostics
	aim[150].gets(aax);
	aim[151].gets(alx);
	aim[152].gets(anx);
}
///////////////////////////////////////////////////////////////////////////////
//Definition of 'intercept' module-variables
//Member function of class 'Aim'
//Module-variable locations are assigned to aim[160-174]
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::def_intercept()
{
	//Definition of module-variables
	aim[161].init("aspazx",0,"Aspect azimuth of incoming missile - deg","intercept","diag","");
	aim[162].init("aspelx",0,"Aspect elevation of incoming missile - deg","intercept","diag","");
}
///////////////////////////////////////////////////////////////////////////////
//'intercept' module 
//Member function of class 'Aim' 
//Parameter Input: 'vehicle_slot' is current 'Aim' object
//Input from module-variable array: 'acft_com_slot' aircraft being attacked, determined in 'seeker' module
//
//Sign convention for incoming missile aspect angles wrt aircraft velocity vector
// Azimuth: positive if missile comes in from the right; negative if from left
// Elevation: positive if missile comes from above, negative if from below;
// Zero is at the positive  direction of the aircraft velocity vector
//
//070412 Created by Peter H Zipfel
///////////////////////////////////////////////////////////////////////////////
void Aim::intercept(Packet *combus,int vehicle_slot,double int_step,char *title)
{
	//local module-variables
	Matrix TTL(3,3);//T.M. of aircraft wrt local level coordinates

	//localizing module-variables
	double aspazx(0);
	double aspelx(0);

	//input data
	//input from other modules
	double time=flat3[0].real();
	Matrix VAEL=flat3[27].vec();//missile
	int acft_com_slot=aim[1].integer();
	Matrix VTEL=aim[2].vec();//aircraft
	double psivlx_acft=aim[3].real();
	double thtvlx_acft=aim[4].real();
	double dta=aim[80].real();
	double dvta=aim[81].real();
	Matrix STAL=aim[89].vec();
	//-------------------------------------------------------------------------
	// displaying miss distance only if missile is inside sphere of aircraft
	if(dta<100){
		// point of closest approach
		if(dvta>0){
			//calculating aspect angles of incoming missile
			//aircraft T velocity relative to incoming missile A
			Matrix VTAEL=VTEL-VAEL;
			//differential speed of missile wrt aircraft
			double diff_speed=VTAEL.absolute();
			//T.M. of aircraft wrt local level coordinates
			TTL=mat2tr(psivlx_acft*RAD,thtvlx_acft*RAD);
			Matrix VTAT=TTL*VTAEL;
			Matrix POLAR=VTAT.pol_from_cart();
			//aspect angles
			double az=POLAR.get_loc(1,0);
			double el=POLAR.get_loc(2,0);
			aspazx=az*DEG;
			aspelx=el*DEG;

			//getting missile # and aircraft #
			string id_acft=combus[acft_com_slot].get_id();
			string id_aim=combus[vehicle_slot].get_id();

			//missile intercepted aircraft
			cout<<"\n"<<" $$$ Intercept of Missile_"<<id_aim<<" with Aircraft_"<<id_acft
				<<"   at sim_time = "<<time<<" sec $$$\n";
			cout<<"      miss distance = "<<dta<<" m     differential speed = "<<diff_speed<<" m/s \n";
			cout<<"      incoming missile azimuth = "<<aspazx<<" deg          elevation = "<<aspelx<<" deg \n\n";

			//missile and aircraft are set to be 'dead'
			combus[vehicle_slot].set_status(0);
			combus[acft_com_slot].set_status(0);
		}
	}	
}
