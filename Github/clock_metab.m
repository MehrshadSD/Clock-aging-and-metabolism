
% read experimental data (used in parameter fitting)

% Specify which simulation to do
% simulation = 1;  % fed (insulin=1, AA=1), no perturbation
% simulation = 2;  % fasting (insulin=0.1, AA=0.5), no perturbation

bmal1=readtable('bmal1.txt');
bmal1.Properties.VariableNames{'Var1'} = 'time';
bmal1.Properties.VariableNames{'Var2'} = 'value';
bmal1=table2array(bmal1);
bmalc = spline(bmal1(:,1),bmal1(:,2),[1:48])';

cry_data=readtable('cry.txt');
cry_data.Properties.VariableNames{'Var1'} = 'time';
cry_data.Properties.VariableNames{'Var2'} = 'value';
cry_data=table2array(cry_data);
cry_data2 = spline(cry_data(:,1),cry_data(:,2),[1:48])';

per_data=readtable('per.txt');
per_data.Properties.VariableNames{'Var1'} = 'time';
per_data.Properties.VariableNames{'Var2'} = 'value';
per_data=table2array(per_data);
per_data2 = spline(per_data(:,1),per_data(:,2),[1:48])';

ror_data=readtable('ROR.txt');
ror_data.Properties.VariableNames{'Var1'} = 'time';
ror_data.Properties.VariableNames{'Var2'} = 'value';
ror_data=table2array(ror_data);
ror_data2 = spline(ror_data(:,1),ror_data(:,2),[1:48])';

rev_data=readtable('REV.txt');
rev_data.Properties.VariableNames{'Var1'} = 'time';
rev_data.Properties.VariableNames{'Var2'} = 'value';
rev_data=table2array(rev_data);
rev_data2 = spline(rev_data(:,1),rev_data(:,2),[1:48])';
 
    
mtorc2_scaling =1;
mtorc1_scaling =1;
glucose_tolerance =0;
base_fig_num =0;
b_pras_mtorc1 =1;
akt_scaling =1;
pi3K_pdk_scaling = 1;
% K_BinhibitM=0.03;
% k_CBM=0.0242;
Nam_transporter = 100;
% K_PM=0.01;  
K_BinhibitM=0.015599155018510 * 5;
K_PM=0.019572951609150;

    load('init_base');  

PP = [0.0104    0.5404    0.2252    0.0151    0.0066    0.0261    0.0849 ...
    0.0390    0.0118   13.2884    0.0261    0.6536   56.1535   14.0871 ...
    3.9968    0.0641   0.3677    0.3005    0.8032];

	tspan=0:1:500;
	opts = odeset('AbsTol',1e-3);
	[t,x]=ode23s(@(t,x) f(t,x,b_pras_mtorc1,glucose_tolerance,mtorc2_scaling,PP,K_BinhibitM,K_PM,simulation,aging,stac_params,nad_params),tspan,x0,opts);
    
% Depending on whether you are using Octave or Matlab,
% you should comment / uncomment one of the following blocks.
% This should also be done for the definition of the function f below.
% Start Matlab code
function xdot=f(t,x,b_pras_mtorc1,glucose_tolerance,mtorc2_scaling,PP,K_BinhibitM,K_PM,simulation,aging,stac_params,nad_params)

    K_bm=PP(1);
    kp_bmal=PP(2);
    dp_bmal=PP(3);
    kass_cb=PP(4);
    kdiss_cb=PP(5);
    k_CBM=PP(6);
    d_cb=PP(7); 
% K_BinhibitM=PP(8);
% K_PM=PP(9);

%     kass_pc=PP(10);
%     kdiss_pc=PP(11);
%     dp_cry=PP(12);
%     dp_per=PP(13);
%     kp_cry=PP(14);
%     kp_per=PP(15); 

% Compartment: id = Cell, name = Cell, constant
	compartment_Cell=1.0;
% Parameter:   id =  IRS_phos_by_Amino_Acids, name = IRS_phos_by_Amino_Acids
	global_par_IRS_phos_by_Amino_Acids=0.0331672;
% Parameter:   id =  AMPK_T172_phos_by_Amino_Acids, name = AMPK_T172_phos_by_Amino_Acids
	global_par_AMPK_T172_phos_by_Amino_Acids=17.6284;
% Parameter:   id =  mTORC2_S2481_phos_by_Amino_Acids, name = mTORC2_S2481_phos_by_Amino_Acids
	global_par_mTORC2_S2481_phos_by_Amino_Acids=0.0268658;
% Parameter:   id =  IR_beta_phos_by_Insulin, name = IR_beta_phos_by_Insulin
	global_par_IR_beta_phos_by_Insulin=0.0203796;
% Parameter:   id =  IR_beta_pY1146_dephos, name = IR_beta_pY1146_dephos
	global_par_IR_beta_pY1146_dephos=0.493514;
% Parameter:   id =  IR_beta_ready, name = IR_beta_ready
	global_par_IR_beta_ready=323.611;
% Parameter:   id =  IRS_phos_by_IR_beta_pY1146, name = IRS_phos_by_IR_beta_pY1146
	global_par_IRS_phos_by_IR_beta_pY1146=2.11894;
% Parameter:   id =  IRS_p_phos_by_p70_S6K_pT229_pT389, name = IRS_p_phos_by_p70_S6K_pT229_pT389
	global_par_IRS_p_phos_by_p70_S6K_pT229_pT389=0.338859859949792;
% Parameter:   id =  IRS_phos_by_p70_S6K_pT229_pT389, name = IRS_phos_by_p70_S6K_pT229_pT389
	global_par_IRS_phos_by_p70_S6K_pT229_pT389=0.0863775267376444;
% Parameter:   id =  IRS_pS636_turnover, name = IRS_pS636_turnover
	global_par_IRS_pS636_turnover=25.0;
% Parameter:   id =  AMPK_T172_phos, name = AMPK_T172_phos
	global_par_AMPK_T172_phos=0.490602;
% Parameter:   id =  AMPK_pT172_dephos, name = AMPK_pT172_dephos
	global_par_AMPK_pT172_dephos=165.704;
% Parameter:   id =  Akt_S473_phos_by_mTORC2_pS2481_first, name = Akt_S473_phos_by_mTORC2_pS2481_first
	global_par_Akt_S473_phos_by_mTORC2_pS2481_first=1.31992E-5;
% Parameter:   id =  Akt_S473_phos_by_mTORC2_pS2481_second, name = Akt_S473_phos_by_mTORC2_pS2481_second
	global_par_Akt_S473_phos_by_mTORC2_pS2481_second=0.159093;
% Parameter:   id =  Akt_T308_phos_by_PI3K_p_PDK1_first, name = Akt_T308_phos_by_PI3K_p_PDK1_first
	global_par_Akt_T308_phos_by_PI3K_p_PDK1_first=7.47437;
% Parameter:   id =  Akt_T308_phos_by_PI3K_p_PDK1_second, name = Akt_T308_phos_by_PI3K_p_PDK1_second
	global_par_Akt_T308_phos_by_PI3K_p_PDK1_second=7.47345;
% Parameter:   id =  Akt_pT308_dephos_first, name = Akt_pT308_dephos_first
	global_par_Akt_pT308_dephos_first=88.9654;
% Parameter:   id =  Akt_pT308_dephos_second, name = Akt_pT308_dephos_second
	global_par_Akt_pT308_dephos_second=88.9639;
% Parameter:   id =  Akt_pS473_dephos_first, name = Akt_pS473_dephos_first
	global_par_Akt_pS473_dephos_first=0.376999;
% Parameter:   id =  Akt_pS473_dephos_second, name = Akt_pS473_dephos_second
	global_par_Akt_pS473_dephos_second=0.380005;
% Parameter:   id =  TSC1_TSC2_S1387_phos_by_AMPK_pT172, name = TSC1_TSC2_S1387_phos_by_AMPK_pT172
	global_par_TSC1_TSC2_S1387_phos_by_AMPK_pT172=0.00175772;
% Parameter:   id =  TSC1_TSC2_T1462_phos_by_Akt_pT308, name = TSC1_TSC2_T1462_phos_by_Akt_pT308
	global_par_TSC1_TSC2_T1462_phos_by_Akt_pT308=1.52417;
% Parameter:   id =  TSC1_TSC2_pS1387_dephos, name = TSC1_TSC2_pS1387_dephos
	global_par_TSC1_TSC2_pS1387_dephos=0.25319;
% Parameter:   id =  TSC1_TSC2_pT1462_dephos, name = TSC1_TSC2_pT1462_dephos
	global_par_TSC1_TSC2_pT1462_dephos=147.239;
% Parameter:   id =  mTORC1_pS2448_dephos_by_TSC1_TSC2, name = mTORC1_pS2448_dephos_by_TSC1_TSC2
	global_par_mTORC1_pS2448_dephos_by_TSC1_TSC2=0.00869774;
% Parameter:   id =  mTORC1_S2448_activation_by_Amino_Acids, name = mTORC1_S2448_activation_by_Amino_Acids
	global_par_mTORC1_S2448_activation_by_Amino_Acids=0.0156992 * 1.5;
% Parameter:   id =  mTORC2_pS2481_dephos, name = mTORC2_pS2481_dephos
	global_par_mTORC2_pS2481_dephos=1.42511;
% Parameter:   id =  mTORC2_S2481_phos_by_PI3K_variant_p, name = mTORC2_S2481_phos_by_PI3K_variant_p
	global_par_mTORC2_S2481_phos_by_PI3K_variant_p=0.120736;
% Parameter:   id =  p70_S6K_T229_phos_by_PI3K_p_PDK1_first, name = p70_S6K_T229_phos_by_PI3K_p_PDK1_first
	global_par_p70_S6K_T229_phos_by_PI3K_p_PDK1_first=0.0133520172873009;
% Parameter:   id =  p70_S6K_T229_phos_by_PI3K_p_PDK1_second, name = p70_S6K_T229_phos_by_PI3K_p_PDK1_second
	global_par_p70_S6K_T229_phos_by_PI3K_p_PDK1_second=1.00000002814509E-6;
% Parameter:   id =  p70_S6K_T389_phos_by_mTORC1_pS2448_first, name = p70_S6K_T389_phos_by_mTORC1_pS2448_first
	global_par_p70_S6K_T389_phos_by_mTORC1_pS2448_first=0.00261303413778722;
% Parameter:   id =  p70_S6K_T389_phos_by_mTORC1_pS2448_second, name = p70_S6K_T389_phos_by_mTORC1_pS2448_second
	global_par_p70_S6K_T389_phos_by_mTORC1_pS2448_second=0.110720890919343;
% Parameter:   id =  p70_S6K_pT229_dephos_first, name = p70_S6K_pT229_dephos_first
	global_par_p70_S6K_pT229_dephos_first=1.00000012897033E-6;
% Parameter:   id =  p70_S6K_pT229_dephos_second, name = p70_S6K_pT229_dephos_second
	global_par_p70_S6K_pT229_dephos_second=0.159201353240651;
% Parameter:   id =  p70_S6K_pT389_dephos_first, name = p70_S6K_pT389_dephos_first
	global_par_p70_S6K_pT389_dephos_first=1.10036057608758;
% Parameter:   id =  p70_S6K_pT389_dephos_second, name = p70_S6K_pT389_dephos_second
	global_par_p70_S6K_pT389_dephos_second=1.10215267954479;
% Parameter:   id =  PRAS40_S183_phos_by_mTORC1_pS2448_first, name = PRAS40_S183_phos_by_mTORC1_pS2448_first
	global_par_PRAS40_S183_phos_by_mTORC1_pS2448_first=0.15881;
% Parameter:   id =  PRAS40_S183_phos_by_mTORC1_pS2448_second, name = PRAS40_S183_phos_by_mTORC1_pS2448_second
	global_par_PRAS40_S183_phos_by_mTORC1_pS2448_second=0.0683009;
% Parameter:   id =  PRAS40_T246_phos_by_Akt_pT308_first, name = PRAS40_T246_phos_by_Akt_pT308_first
	global_par_PRAS40_T246_phos_by_Akt_pT308_first=0.279344;
% Parameter:   id =  PRAS40_T246_phos_by_Akt_pT308_second, name = PRAS40_T246_phos_by_Akt_pT308_second
	global_par_PRAS40_T246_phos_by_Akt_pT308_second=0.279401;
% Parameter:   id =  PRAS40_pS183_dephos_first, name = PRAS40_pS183_dephos_first
	global_par_PRAS40_pS183_dephos_first=1.8706;
% Parameter:   id =  PRAS40_pS183_dephos_second, name = PRAS40_pS183_dephos_second
	global_par_PRAS40_pS183_dephos_second=1.88453;
% Parameter:   id =  PRAS40_pT246_dephos_first, name = PRAS40_pT246_dephos_first
	global_par_PRAS40_pT246_dephos_first=11.8759;
% Parameter:   id =  PRAS40_pT246_dephos_second, name = PRAS40_pT246_dephos_second
	global_par_PRAS40_pT246_dephos_second=11.876;
% Parameter:   id =  PI3K_p_PDK1_dephos, name = PI3K_p_PDK1_dephos
	global_par_PI3K_p_PDK1_dephos=0.18913343080532;
% Parameter:   id =  PI3K_PDK1_phos_by_IRS_p, name = PI3K_PDK1_phos_by_IRS_p
	global_par_PI3K_PDK1_phos_by_IRS_p=1.87226757782201E-4;
% Parameter:   id =  PI3K_variant_p_dephos, name = PI3K_variant_p_dephos
	global_par_PI3K_variant_p_dephos=0.108074886441184;
% Parameter:   id =  PI3K_variant_phos_by_IR_beta_pY1146, name = PI3K_variant_phos_by_IR_beta_pY1146
	global_par_PI3K_variant_phos_by_IR_beta_pY1146=5.49027801822575E-4;
% Parameter:   id =  scale_IR_beta_pY1146_obs, name = scale_IR_beta_pY1146_obs
	global_par_scale_IR_beta_pY1146_obs=1.0;
% Parameter:   id =  scale_IRS_pS636_obs, name = scale_IRS_pS636_obs
	global_par_scale_IRS_pS636_obs=1.0;
% Parameter:   id =  scale_AMPK_pT172_obs, name = scale_AMPK_pT172_obs
	global_par_scale_AMPK_pT172_obs=1.0;
% Parameter:   id =  scale_Akt_pT308_obs, name = scale_Akt_pT308_obs
	global_par_scale_Akt_pT308_obs=1.0;
% Parameter:   id =  scale_Akt_pS473_obs, name = scale_Akt_pS473_obs
	global_par_scale_Akt_pS473_obs=1.0;
% Parameter:   id =  scale_TSC1_TSC2_pS1387_obs, name = scale_TSC1_TSC2_pS1387_obs
	global_par_scale_TSC1_TSC2_pS1387_obs=1.0;
% Parameter:   id =  scale_mTOR_pS2448_obs, name = scale_mTOR_pS2448_obs
	global_par_scale_mTOR_pS2448_obs=1.0;
% Parameter:   id =  scale_mTOR_pS2481_obs, name = scale_mTOR_pS2481_obs
	global_par_scale_mTOR_pS2481_obs=1.0;
% Parameter:   id =  scale_p70_S6K_pT229_obs, name = scale_p70_S6K_pT229_obs
	global_par_scale_p70_S6K_pT229_obs=1.0;
% Parameter:   id =  scale_p70_S6K_pT389_obs, name = scale_p70_S6K_pT389_obs
	global_par_scale_p70_S6K_pT389_obs=1.0;
% Parameter:   id =  scale_PRAS40_pT246_obs, name = scale_PRAS40_pT246_obs
	global_par_scale_PRAS40_pT246_obs=1.0;
% Parameter:   id =  scale_PRAS40_pS183_obs, name = scale_PRAS40_pS183_obs
	global_par_scale_PRAS40_pS183_obs=1.0;
% assignmentRule: variable = Insulin
% 	x(32)=piecewise(0, time < 0, 0);
    if simulation == 1  % fed
        x(32)=1;
        x(33)=1;
    elseif simulation == 2  % fasting
        x(32)=0.1;
        x(33)=0.5;
    elseif simulation == 3  % eat half the day
%         if t <= 117+12 | (t > 117+24 & t <= 117+36)
%             x(32)=1;
%             x(33)=1;
%         else
%             x(32)=0.1;
%             x(33)=0.5;
%         end
        x(32) = 0.1+0.9*get_P(t,6,12,0.1,24);
        x(33) = 0.5+0.5*get_P(t,6,12,0.1,24);
    elseif simulation == 4  % eat the other half the day
%         if t >= 117+12
%             x(32)=1;
%             x(33)=1;
%         else
%             x(32)=0.1;
%             x(33)=0.5;
%         end
        x(32) = 0.1+0.9*get_P(t,18,12,0.1,24);
        x(33) = 0.5+0.5*get_P(t,18,12,0.1,24);
%         x(32) = 0.1+0.9*get_P(t,16,12,0.1,24);
%         x(33) = 0.5+0.5*get_P(t,16,12,0.1,24);
    else
        display('simulation not known');
    end
% assignmentRule: variable = IR_beta_pY1146_obs
	x(34)=x(2);
% assignmentRule: variable = IRS_pS636_obs
	x(35)=x(6);
% assignmentRule: variable = AMPK_pT172_obs
	x(36)=x(8);
% assignmentRule: variable = Akt_pT308_obs
	x(37)=x(10)+x(12);
% assignmentRule: variable = Akt_pS473_obs
	x(38)=x(11)+x(12);
% assignmentRule: variable = TSC1_TSC2_pS1387_obs
	x(39)=x(15);
% assignmentRule: variable = mTOR_pS2448_obs
	x(40)=x(17);
% assignmentRule: variable = mTOR_pS2481_obs
	x(41)=x(19);
% assignmentRule: variable = p70_S6K_pT229_obs
	x(42)=x(21)+x(23);
% assignmentRule: variable = p70_S6K_pT389_obs
	x(43)=x(22)+x(23);
% assignmentRule: variable = PRAS40_pT246_obs
	x(44)=x(25)+x(27);
% assignmentRule: variable = PRAS40_pS183_obs
	x(45)=x(26)+x(27);

    
    kUKAP= 0.000163;

    ksirt=5*2.0e4;
%     ksirt=5;


% Retrieve variables
%     NAMPT = x(46);
%     NAM   = x(47);
%     NMN   = x(48);
%     NADp  = x(49);
    SIRT  = x(56);  % made obsolute because now we use values computed in the clock
    
    FOXO  = x(57);
    FOXOc = x(58);
    FOXOp = x(59);
    PGC1  = x(60);
    ULK1s = x(61);
    ULK   = x(62);
    
    oldNAMPT = x(63);
    oldNAM   = x(64);
    oldNMN   = x(65);
    oldNADp  = x(66);
    oldSIRT  = x(67);
    
    xdot=zeros(84,1);
    
    ampk_scaling = 1.2/28.35; % LeFranc fast+fed average / our AMPK fast+fed average
    clock_vars = circadian_mtorc(t,x(69:84),x(17),x(8)*ampk_scaling,PP,aging,stac_params,nad_params);
    xdot(69:84) = clock_vars(1:end-1);
    SIRT  = clock_vars(end);

% Reaction: id = reaction_1, name = reaction_1
	reaction_reaction_1=compartment_Cell*function_1(x(33), x(4), global_par_IRS_phos_by_Amino_Acids);

% Reaction: id = reaction_2, name = reaction_2
	reaction_reaction_2=compartment_Cell*function_2(x(7), global_par_AMPK_T172_phos_by_Amino_Acids, x(33));

% Reaction: id = reaction_3, name = reaction_3
	reaction_reaction_3=compartment_Cell*function_3(x(33), x(18), global_par_mTORC2_S2481_phos_by_Amino_Acids);

% Reaction: id = reaction_4, name = reaction_4
	reaction_reaction_4=compartment_Cell*function_4(x(1), global_par_IR_beta_phos_by_Insulin, x(32));

% Reaction: id = reaction_5, name = reaction_5
	reaction_reaction_5=compartment_Cell*global_par_IR_beta_pY1146_dephos*x(2);

% Reaction: id = reaction_6, name = reaction_6
	reaction_reaction_6=compartment_Cell*global_par_IR_beta_ready*x(3);

% Reaction: id = reaction_7, name = reaction_7
	reaction_reaction_7=compartment_Cell*function_5(x(4), global_par_IRS_phos_by_IR_beta_pY1146, x(2));

% Reaction: id = reaction_8, name = reaction_8
	reaction_reaction_8=compartment_Cell*function_6(x(5), global_par_IRS_p_phos_by_p70_S6K_pT229_pT389, x(23));

% Reaction: id = reaction_9, name = reaction_9
	reaction_reaction_9=compartment_Cell*function_7(x(4), global_par_IRS_phos_by_p70_S6K_pT229_pT389, x(23));

% Reaction: id = reaction_10, name = reaction_10
	reaction_reaction_10=compartment_Cell*global_par_IRS_pS636_turnover*x(6);

% Reaction: id = reaction_11, name = reaction_11
	reaction_reaction_11=compartment_Cell*global_par_PI3K_p_PDK1_dephos*x(31);

% Reaction: id = reaction_12, name = reaction_12
	reaction_reaction_12=compartment_Cell*function_8(x(5), x(30), global_par_PI3K_PDK1_phos_by_IRS_p);

% Reaction: id = reaction_13, name = reaction_13
	reaction_reaction_13=compartment_Cell*global_par_PI3K_variant_p_dephos*x(29);

% Reaction: id = reaction_14, name = reaction_14
	reaction_reaction_14=compartment_Cell*function_9(x(2), x(28), global_par_PI3K_variant_phos_by_IR_beta_pY1146);

% Reaction: id = reaction_15, name = reaction_15
	reaction_reaction_15=compartment_Cell*function_10(x(7), global_par_AMPK_T172_phos, x(5));

% Reaction: id = reaction_16, name = reaction_16
	reaction_reaction_16=compartment_Cell*global_par_AMPK_pT172_dephos*x(8);

% Reaction: id = reaction_17, name = reaction_17
	reaction_reaction_17=compartment_Cell*function_11(x(9), global_par_Akt_T308_phos_by_PI3K_p_PDK1_first, x(31));

% Reaction: id = reaction_18, name = reaction_18
	reaction_reaction_18=compartment_Cell*function_12(x(9), global_par_Akt_S473_phos_by_mTORC2_pS2481_first, x(19));

% Reaction: id = reaction_19, name = reaction_19
	reaction_reaction_19=compartment_Cell*function_13(global_par_Akt_T308_phos_by_PI3K_p_PDK1_second, x(11), x(31));

% Reaction: id = reaction_20, name = reaction_20
	reaction_reaction_20=compartment_Cell*function_14(global_par_Akt_S473_phos_by_mTORC2_pS2481_second, x(10), x(19));

% Reaction: id = reaction_21, name = reaction_21
	reaction_reaction_21=compartment_Cell*global_par_Akt_pT308_dephos_first*x(10);

% Reaction: id = reaction_22, name = reaction_22
	reaction_reaction_22=compartment_Cell*global_par_Akt_pS473_dephos_first*x(11);

% Reaction: id = reaction_23, name = reaction_23
	reaction_reaction_23=compartment_Cell*global_par_Akt_pT308_dephos_second*x(12);

% Reaction: id = reaction_24, name = reaction_24
	reaction_reaction_24=compartment_Cell*global_par_Akt_pS473_dephos_second*x(12);

% Reaction: id = reaction_25, name = reaction_25
	reaction_reaction_25=compartment_Cell*function_15(x(8), x(13), global_par_TSC1_TSC2_S1387_phos_by_AMPK_pT172);

% Reaction: id = reaction_26, name = reaction_26
	reaction_reaction_26=compartment_Cell*function_16(x(10), x(12), x(13), global_par_TSC1_TSC2_T1462_phos_by_Akt_pT308);

% Reaction: id = reaction_27, name = reaction_27
	reaction_reaction_27=compartment_Cell*global_par_TSC1_TSC2_pS1387_dephos*x(15);

% Reaction: id = reaction_28, name = reaction_28
	reaction_reaction_28=compartment_Cell*global_par_TSC1_TSC2_pT1462_dephos*x(14);

% Reaction: id = reaction_29, name = reaction_29
	reaction_reaction_29=compartment_Cell*function_17(x(13), x(15), x(17), global_par_mTORC1_pS2448_dephos_by_TSC1_TSC2);

% Reaction: id = reaction_30, name = reaction_30
	reaction_reaction_30=compartment_Cell*function_18(x(33), x(16), global_par_mTORC1_S2448_activation_by_Amino_Acids);

% Reaction: id = reaction_31, name = reaction_31
	reaction_reaction_31=compartment_Cell*function_19(x(29), x(18), global_par_mTORC2_S2481_phos_by_PI3K_variant_p);

% Reaction: id = reaction_32, name = reaction_32
	reaction_reaction_32=compartment_Cell*global_par_mTORC2_pS2481_dephos*x(19);

% Reaction: id = reaction_33, name = reaction_33
	reaction_reaction_33=compartment_Cell*function_20(x(31), x(20), global_par_p70_S6K_T229_phos_by_PI3K_p_PDK1_first);

% Reaction: id = reaction_34, name = reaction_34
	reaction_reaction_34=compartment_Cell*function_21(x(17), x(20), global_par_p70_S6K_T389_phos_by_mTORC1_pS2448_first);

% Reaction: id = reaction_35, name = reaction_35
	reaction_reaction_35=compartment_Cell*function_22(x(31), global_par_p70_S6K_T229_phos_by_PI3K_p_PDK1_second, x(22));

% Reaction: id = reaction_36, name = reaction_36
	reaction_reaction_36=compartment_Cell*function_23(x(17), global_par_p70_S6K_T389_phos_by_mTORC1_pS2448_second, x(21));

% Reaction: id = reaction_37, name = reaction_37
	reaction_reaction_37=compartment_Cell*global_par_p70_S6K_pT229_dephos_first*x(21);

% Reaction: id = reaction_38, name = reaction_38
	reaction_reaction_38=compartment_Cell*global_par_p70_S6K_pT389_dephos_first*x(22);

% Reaction: id = reaction_39, name = reaction_39
	reaction_reaction_39=compartment_Cell*global_par_p70_S6K_pT229_dephos_second*x(23);

% Reaction: id = reaction_40, name = reaction_40
	reaction_reaction_40=compartment_Cell*global_par_p70_S6K_pT389_dephos_second*x(23);

% Reaction: id = reaction_41, name = reaction_41
	reaction_reaction_41=compartment_Cell*function_24(x(24), global_par_PRAS40_S183_phos_by_mTORC1_pS2448_first, x(17));

% Reaction: id = reaction_42, name = reaction_42
	reaction_reaction_42=compartment_Cell*function_25(x(10), x(12), x(24), global_par_PRAS40_T246_phos_by_Akt_pT308_first);

% Reaction: id = reaction_43, name = reaction_43
	reaction_reaction_43=compartment_Cell*function_26(x(10), x(12), global_par_PRAS40_T246_phos_by_Akt_pT308_second, x(26));

% Reaction: id = reaction_44, name = reaction_44
	reaction_reaction_44=compartment_Cell*function_27(global_par_PRAS40_S183_phos_by_mTORC1_pS2448_second, x(25), x(17));

% Reaction: id = reaction_45, name = reaction_45
	reaction_reaction_45=compartment_Cell*global_par_PRAS40_pS183_dephos_first*x(26);

% Reaction: id = reaction_46, name = reaction_46
	reaction_reaction_46=compartment_Cell*global_par_PRAS40_pT246_dephos_first*x(25);

% Reaction: id = reaction_47, name = reaction_47
	reaction_reaction_47=compartment_Cell*global_par_PRAS40_pS183_dephos_second*x(27);

% Reaction: id = reaction_48, name = reaction_48
	reaction_reaction_48=compartment_Cell*global_par_PRAS40_pT246_dephos_second*x(27);
	
% Species:   id = IR_beta, name = IR_beta, affected by kineticLaw
	xdot(1) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_4) + ( 1.0 * reaction_reaction_6));
	
% Species:   id = IR_beta_pY1146, name = IR_beta_pY1146, affected by kineticLaw
	xdot(2) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_4) + (-1.0 * reaction_reaction_5));
	
% Species:   id = IR_beta_refractory, name = IR_beta_refractory, affected by kineticLaw
	xdot(3) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_5) + (-1.0 * reaction_reaction_6));
	
% Species:   id = IRS, name = IRS, affected by kineticLaw
	xdot(4) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_1) + (-1.0 * reaction_reaction_7) + (-1.0 * reaction_reaction_9) + ( 1.0 * reaction_reaction_10));
	
% Species:   id = IRS_p, name = IRS_p, affected by kineticLaw
	xdot(5) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_1) + ( 1.0 * reaction_reaction_7) + (-1.0 * reaction_reaction_8));

% Species:   id = IRS_pS636, name = IRS_pS636, affected by kineticLaw
	xdot(6) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_8) + ( 1.0 * reaction_reaction_9) + (-1.0 * reaction_reaction_10));
	
% Species:   id = AMPK, name = AMPK, affected by kineticLaw
	xdot(7) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_2) + (-1.0 * reaction_reaction_15) + ( 1.0 * reaction_reaction_16))-ksirt*SIRT*x(7) + (ULK1s*x(8)*kUKAP);
	
% Species:   id = AMPK_pT172, name = AMPK_pT172, affected by kineticLaw
	xdot(8) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_2) + ( 1.0 * reaction_reaction_15) + (-1.0 * reaction_reaction_16))+ ksirt*SIRT*x(7) - (ULK1s*x(8)*kUKAP);
	
% [xdot(8),( 1.0 * reaction_reaction_2),( 1.0 * reaction_reaction_15),(-1.0 * reaction_reaction_16), ksirt*SIRT*x(7), - (ULK1s*x(8)*kUKAP), (1/(compartment_Cell))]
% pause
    
% Species:   id = Akt, name = Akt, affected by kineticLaw
	xdot(9) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_17) + (-1.0 * reaction_reaction_18) + ( 1.0 * reaction_reaction_21) + ( 1.0 * reaction_reaction_22));
	
% Species:   id = Akt_pT308, name = Akt_pT308, affected by kineticLaw
	xdot(10) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_17) + (-1.0 * reaction_reaction_20) + (-1.0 * reaction_reaction_21) + ( 1.0 * reaction_reaction_24));
	
% Species:   id = Akt_pS473, name = Akt_pS473, affected by kineticLaw
	xdot(11) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_18) + (-1.0 * reaction_reaction_19) + (-1.0 * reaction_reaction_22) + ( 1.0 * reaction_reaction_23));
	
% Species:   id = Akt_pT308_pS473, name = Akt_pT308_pS473, affected by kineticLaw
	xdot(12) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_19) + ( 1.0 * reaction_reaction_20) + (-1.0 * reaction_reaction_23) + (-1.0 * reaction_reaction_24));
	
% Species:   id = TSC1_TSC2, name = TSC1_TSC2, affected by kineticLaw
	xdot(13) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_25) + (-1.0 * reaction_reaction_26) + ( 1.0 * reaction_reaction_27) + ( 1.0 * reaction_reaction_28));
	
% Species:   id = TSC1_TSC2_pT1462, name = TSC1_TSC2_pT1462, affected by kineticLaw
	xdot(14) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_26) + (-1.0 * reaction_reaction_28));
	
% Species:   id = TSC1_TSC2_pS1387, name = TSC1_TSC2_pS1387, affected by kineticLaw
	xdot(15) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_25) + (-1.0 * reaction_reaction_27));
	
% Species:   id = mTORC1, name = mTORC1, affected by kineticLaw
	xdot(16) = (1/(compartment_Cell))*(( 1 * reaction_reaction_29) + (-1.0 * reaction_reaction_30)) + b_pras_mtorc1*0.00068*x(17)*10^(0.25*x(24)) + K_BinhibitM*x(78)*x(17) + K_PM*x(74)*x(17);
% 	xdot(16) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_29) + (-1.0 * reaction_reaction_30));
% x(32), 1 * reaction_reaction_29, + (-1.0 * reaction_reaction_30),b_pras_mtorc1*0.00068*x(17)*10^(0.25*x(24)), + K_BinhibitM*x(78)*x(17), + K_PM*x(74)*x(17)
% pause
	
% Species:   id = mTORC1_pS2448, name = mTORC1_pS2448, affected by kineticLaw
	xdot(17) = (1/(compartment_Cell))*((-1 * reaction_reaction_29) + ( 1.0 * reaction_reaction_30)) - b_pras_mtorc1*0.00068*x(17)*10^(0.25*x(24)) - K_BinhibitM*x(78)*x(17) - K_PM*x(74)*x(17);
% 	xdot(17) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_29) + ( 1.0 * reaction_reaction_30));
	
% Species:   id = mTORC2, name = mTORC2, affected by kineticLaw
	xdot(18) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_3) + (-1.0 * reaction_reaction_31) + ( 1.0 * reaction_reaction_32));
	
% Species:   id = mTORC2_pS2481, name = mTORC2_pS2481, affected by kineticLaw
	xdot(19) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_3) + ( 1.0 * reaction_reaction_31) + (-1.0 * reaction_reaction_32));
	
% Species:   id = p70_S6K, name = p70_S6K, affected by kineticLaw
	xdot(20) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_33) + (-1.0 * reaction_reaction_34) + ( 1.0 * reaction_reaction_37) + ( 1.0 * reaction_reaction_38));
	
% Species:   id = p70_S6K_pT229, name = p70_S6K_pT229, affected by kineticLaw
	xdot(21) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_33) + (-1.0 * reaction_reaction_36) + (-1.0 * reaction_reaction_37) + ( 1.0 * reaction_reaction_40));
	
% Species:   id = p70_S6K_pT389, name = p70_S6K_pT389, affected by kineticLaw
	xdot(22) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_34) + (-1.0 * reaction_reaction_35) + (-1.0 * reaction_reaction_38) + ( 1.0 * reaction_reaction_39));
	
% Species:   id = p70_S6K_pT229_pT389, name = p70_S6K_pT229_pT389, affected by kineticLaw
	xdot(23) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_35) + ( 1.0 * reaction_reaction_36) + (-1.0 * reaction_reaction_39) + (-1.0 * reaction_reaction_40));
	
% Species:   id = PRAS40, name = PRAS40, affected by kineticLaw
	xdot(24) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_41) + (-1.0 * reaction_reaction_42) + ( 1.0 * reaction_reaction_45) + ( 1.0 * reaction_reaction_46));
	
% Species:   id = PRAS40_pT246, name = PRAS40_pT246, affected by kineticLaw
	xdot(25) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_42) + (-1.0 * reaction_reaction_44) + (-1.0 * reaction_reaction_46) + ( 1.0 * reaction_reaction_47));
	
% Species:   id = PRAS40_pS183, name = PRAS40_pS183, affected by kineticLaw
	xdot(26) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_41) + (-1.0 * reaction_reaction_43) + (-1.0 * reaction_reaction_45) + ( 1.0 * reaction_reaction_48));
	
% Species:   id = PRAS40_pT246_pS183, name = PRAS40_pT246_pS183, affected by kineticLaw
	xdot(27) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_43) + ( 1.0 * reaction_reaction_44) + (-1.0 * reaction_reaction_47) + (-1.0 * reaction_reaction_48));
	
% Species:   id = PI3K_variant, name = PI3K_variant, affected by kineticLaw
	xdot(28) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_13) + (-1.0 * reaction_reaction_14));
	
% Species:   id = PI3K_variant_p, name = PI3K_variant_p, affected by kineticLaw
	xdot(29) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_13) + ( 1.0 * reaction_reaction_14));
	
% Species:   id = PI3K_PDK1, name = PI3K_PDK1, affected by kineticLaw
	xdot(30) = (1/(compartment_Cell))*(( 1.0 * reaction_reaction_11) + (-1.0 * reaction_reaction_12));
	
% Species:   id = PI3K_p_PDK1, name = PI3K_p_PDK1, affected by kineticLaw
	xdot(31) = (1/(compartment_Cell))*((-1.0 * reaction_reaction_11) + ( 1.0 * reaction_reaction_12));
	
    
    
%LIVER CLOCK LEFRANC  Kp_NAMPT*(NAMPT)=70
%NAMPT
m_NAMPT_AMPK=0.6352243791;
dp_NAMPT=49.8841023982;
km4=20;
kkm4=20;
km5=40;
% km5=4;
kkm5=5;
km2=40;
kkm2=1;
km1=5;
kkm1=2;
km6=5;
kmm6=1;
km7=2;
kmm7=1;
NAD=10;
% SIRTT=5;

%     NAM   = x(53);
%     NMN   = x(50);
%     NADp  = x(54);
%     SIRT  = x(56);
    
%Protein NAMPT
%     xdot(46)=50 - ((dp_NAMPT*NAMPT)/(1+m_NAMPT_AMPK*x(8)));
% 
% % circadian clock
% %Protein NAM
%     xdot(53)=(km4*NADp/(kkm4+NADp)) - ((km5*NAMPT*NAM)/(kkm5+NAM));
% 
% %Protein NMN
%     xdot(50)=(km5*NAMPT*NAM/(kkm5+NAM)) - ((km2*NMN)/(kkm2+NMN));
% 
% %protein NADP
%     xdot(54)=km1*(NAD-NADp)/(kkm1+(NAD-NADp)) + (km2*NMN/(kkm2+NMN)) - (km4*(NADp)/kkm4*NADp);
% 
% %SIRT
%     sirt0=(150*NADp)./(251+NADp);
%     xdot(56) = 1000*(sirt0-SIRT);
    
    xdot(46:56) = NADSOLV(t*60,x(46:56),x(8),x(83));   % converted time from min to sec
    xdot(56) = 100*(SIRT-x(56));
   
    
%old NAMPT
    xdot(63)=50 - ((dp_NAMPT*oldNAMPT)/(1+m_NAMPT_AMPK*x(8)));
%old NAM
    xdot(64)=(km4*oldNADp/(kkm4+oldNADp)) - ((km5*oldNAMPT*oldNAM)/(kkm5+oldNAM));
%oid NMN
    xdot(65)=(km5*oldNAMPT*oldNAM/(kkm5+oldNAM)) - ((km2*oldNMN)/(kkm2+oldNMN));
%protein NADP
    xdot(66)=km1*(NAD-oldNADp)/(kkm1+(NAD-oldNADp)) + (km2*oldNMN/(kkm2+oldNMN)) - (km4*(oldNADp)/(kkm4+oldNADp));
%SIRT
    sirt0=(150*oldNADp)./(251+oldNADp);
    xdot(67) = 1000*(sirt0-oldSIRT);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLASMA GLUCOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if (240<t & t<240+25)
% %         kglu_in    = 1.2e-1*(25-t)/25;
%         kglu_in    = 1e-1;
%     else
%         kglu_in    = 0;
%     end
%     kglu_metab = 0.0001;
%     xdot(68) = kglu_in - kglu_metab*x(68)*x(12)^2;
    
%     if (240<t & t<240+25)
% %         kglu_in    = 1.2e-1*(25-t)/25;
%         kglu_in    = 0.115-(t-240)/25*0.115;
%     else
%         kglu_in    = 0;
%     end
    kglu_in    = 0;
    kglu_metab = 0.0001*7;
    xdot(68) = kglu_in - kglu_metab*x(68)*x(12);
    
%     ampk_scaling = 1.2/28.35; % LeFranc fast+fed average / our AMPK fast+fed average
%     xdot(69:84) = circadian_mtorc(t,x(69:84),x(17),x(8)*ampk_scaling,PP);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FOXO DYNAMIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%

kAKT=0.5; %5*10^-5*10^4 

kAMPK=0.5;   %5*10^-5*10^4 

kSIRTFOXO=0.1;  %10^-4*10^3 

kCBp=10^-4;  %10^-4*10^3 

CBP=10^3;

    xdot(57)=kAKT*(x(10)+x(12))*(FOXOc) - kAMPK*(x(8))*(FOXO);  %FOXO
   
    xdot(58)=kAMPK*(FOXO)*(x(8)) - kAKT*(x(12)+x(10))*(FOXOc) - kSIRTFOXO*(FOXOc)*(oldSIRT) + kCBp*(FOXOp)*(CBP);  %FOXO_C

    xdot(59)=kSIRTFOXO*(FOXOc)*(oldSIRT) - kCBp*(FOXOp)*(CBP); %FOXO_P

%%%%%%%%%%%%%%%%%%%%%%%%%%% PGC1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
Vopgc=0.046;
Vspgc=0.0142;
kAPGC=7.578;
Vdpgc=0.388;
kdpgc=3.299;
kdn=0.077;

    xdot(60)=Vopgc + (Vspgc*(oldSIRT)/(kAPGC+oldSIRT)) - (Vdpgc*PGC1/(kdpgc+PGC1)) - kdn*PGC1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ULK1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k10 = 0.000166;
kULKd = 0.00169;
KULKM = 0.000158;

%ULk1S
    xdot(61)= k10*ULK*x(8) - kULKd*ULK1s - KULKM*ULK1s*x(17);

%ULK
    xdot(62)= -k10*ULK*x(8) + kULKd*ULK1s + KULKM*ULK1s*x(17);
 
   
% Species:   id = Insulin, name = Insulin, involved in a rule 	xdot(32) = x(32);
	
% Species:   id = Amino_Acids, name = Amino_Acids, involved in a rule 	xdot(33) = x(33);
	
% Species:   id = IR_beta_pY1146_obs, name = IR_beta_pY1146_obs, involved in a rule 	xdot(34) = x(34);
	
% Species:   id = IRS_pS636_obs, name = IRS_pS636_obs, involved in a rule 	xdot(35) = x(35);
	
% Species:   id = AMPK_pT172_obs, name = AMPK_pT172_obs, involved in a rule 	xdot(36) = x(36);
	
% Species:   id = Akt_pT308_obs, name = Akt_pT308_obs, involved in a rule 	xdot(37) = x(37);
	
% Species:   id = Akt_pS473_obs, name = Akt_pS473_obs, involved in a rule 	xdot(38) = x(38);
	
% Species:   id = TSC1_TSC2_pS1387_obs, name = TSC1_TSC2_pS1387_obs, involved in a rule 	xdot(39) = x(39);
	
% Species:   id = mTOR_pS2448_obs, name = mTOR_pS2448_obs, involved in a rule 	xdot(40) = x(40);
	
% Species:   id = mTOR_pS2481_obs, name = mTOR_pS2481_obs, involved in a rule 	xdot(41) = x(41);
	
% Species:   id = p70_S6K_pT229_obs, name = p70_S6K_pT229_obs, involved in a rule 	xdot(42) = x(42);
	
% Species:   id = p70_S6K_pT389_obs, name = p70_S6K_pT389_obs, involved in a rule 	xdot(43) = x(43);
	
% Species:   id = PRAS40_pT246_obs, name = PRAS40_pT246_obs, involved in a rule 	xdot(44) = x(44);
	
% Species:   id = PRAS40_pS183_obs, name = PRAS40_pS183_obs, involved in a rule 	xdot(45) = x(45);

%save('mainx.mat','xdot'), pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Function for computing the N-species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T =NADSOLV(t,x,pAMPK,NADT)

%   pAMPK will increase NAMPT, which eventually increases SIRT=

     Nam_ex = x(1);
     NADbound = x(2);
     NR = x(3);
     NaAD = x(4);
     NMN = x(5);
     NA = x(6);
     NaMN = x(7);
     NAM = x(8);
     NAD = x(9);
     NAR = x(10);
     mySIRT = x(11);

%      Nam_ex = x(46);
%      NADbound = x(47);
%      NR = x(48);
%      NaAD = x(49);
%      NMN = x(50);
%      NA = x(51);
%      NaMN = x(52);
%      NAM = x(53);
%      NAD = x(54);
%      NAR = x(55);
%      mySIRT = x(56);
  

 
  
     
%%   Parameter Values:
     Nam_transporter = 100;
     Km_NamPRT = 5e-6;
     Kcat_NamPRT = 0.0077;
%      cell_devision_all = 2.8e-05;
     ETNMNAT = 10;
     Initial_for_cell_devision_all = 2.8e-05;
     NAPRT_E_T = 10;
     NAPRT_turnover = 3.3;
     NAPRT_scaling = 0.001;
     NAPRT_Km = 0.00015;
     NMNAT1_NMN_scaling = 0.001;
     NMNAT1_NMN_kcat_A = 53.8;
     NMNAT1_NMN_Km_A = 0.0223;
     NMNAT1_NMN_kcat_PA = 129.1;
     NMNAT1_NMN_Km_PA = 0.059;
     NMNAT1_NMN_Km_B = 0.0677;
     NMNAT1_NMN_Km_PB = 0.502;
     NMNAT1_NaMN_scaling = 0.001;
     NMNAT1_NaMN_kcat_A = 42.9;
     NMNAT1_NaMN_Km_A = 0.0677;
     NMNAT1_NaMN_kcat_PA = 103.8;
     NMNAT1_NaMN_Km_PA = 0.502;
     NMNAT1_NaMN_Km_B = 22.3;
     NMNAT1_NaMN_Km_PB = 0.059;
     NADS_E_T = 10;
     NADS_turnover = 21;
     NADS_scaling = 0.001;
     NADS_Km = 0.19;
     NT5_NMN_ET = 10;
     NT5_NMN_scaling = 0.001;
     NT5_NMN_kcat = 0.5;
     NT5_NMN_Km = 5;
     PNP_NR_E_T = 10;
     PNP_NR_turnover = 40;
     PNP_NR_scaling = 0.001;
     PNP_NR_Km = 1.48;
     NNMT_ET = 10;
     NNMT_scaling = 0.001;
     NNMT_Kcat = 8.1;
     NNMT_Ki = 0.06;
     Kma = 0.0018;
     Kmb = 0.4;
     NT5_NaMN_ET = 10;
     NT5_NaMN_scaling = 0.001;
     NT5_NaMN_kcat = 2.8;
     NT5_NaMN_Km = 3.5;
     SIRT_ET = 10;
     SIRT_scaling = 0.001;
     SIRT_Kcat = 0.67;
     SIRT_Km = 0.029;
     SIRT_Ki = 0.06;
     NAD_consumption_without_Nam_inhibition_ET = 10;
     NAD_consumption_without_Nam_inhibition_scaling = 0.001;
     NAD_consumption_without_Nam_inhibition_kcat = 1;
     NAD_consumption_without_Nam_inhibition_Km = 1;
     NADA_ET = 0;
     NADA_scaling = 0.001;
     NADA_Kcat = 0.69;
     NADA_Km = 0.009;
     NADA_Ki = 0.1;
     NAMPT_ET = 400;
     NAMPT_scaling = 0.001;
     NAMPT_Ki = 0.0021;
     NAD_binding_k1 = 100;
     NAD_binding_k2 = 10;
     NRK1_NaMN_ET = 10;
     NRK1_NaMN_scaling = 0.001;
     NRK1_NaMN_kcat = 0.23;
     NRK1_NaMN_Km = 0.0034;
     PNP_NAR_E_T = 10;
     PNP_NAR_turnover = 40;
     PNP_NAR_scaling = 0.001;
     PNP_NAR_Km = 1.48;
     NRK1_NMN_ET = 10;
     NRK1_NMN_scaling = 0.001;
     NRK1_NMN_kcat = 0.23;
     NRK1_NMN_Km = 0.0034;


     
     
     
     %      NAPRT(1)).E_T = 10
%      NAPRT(1)).turnover = 3.3
%      NAPRT(1)).scaling = 0.001
%      NAPRT(1)).Km = 0.00015
%      NMNAT1_NMN(1)).scaling = 0.001
%      NMNAT1_NMN(1)).kcat_A = 53.8
%      NMNAT1_NMN(1)).Km_A = 0.0223
%      (NMNAT1_NMN(1)).kcat_PA = 129.1
%      (NMNAT1_NMN(1)).Km_PA = 0.059
%      (NMNAT1_NMN(1)).Km_B = 0.0677
%      (NMNAT1_NMN(1)).Km_PB = 0.502
%      (NMNAT1_NaMN(1)).scaling = 0.001
%      (NMNAT1_NaMN(1)).kcat_A = 42.9
%      (NMNAT1_NaMN(1)).Km_A = 0.0677
%      (NMNAT1_NaMN(1)).kcat_PA = 103.8
%      (NMNAT1_NaMN(1)).Km_PA = 0.502
%      (NMNAT1_NaMN(1)).Km_B = 22.3
%      (NMNAT1_NaMN(1)).Km_PB = 0.059
%      (NADS(1)).E_T = 10
%      (NADS(1)).turnover = 21
%      (NADS(1)).scaling = 0.001
%      (NADS(1)).Km = 0.19
%      (NT5_NMN(1)).ET = 10
%      (NT5_NMN(1)).scaling = 0.001
%      (NT5_NMN(1)).kcat = 0.5
%      (NT5_NMN(1)).Km = 5
%      (PNP_NR(1)).E_T = 10
%      (PNP_NR(1)).turnover = 40
%      (PNP_NR(1)).scaling = 0.001
%      (PNP_NR(1)).Km = 1.48
%      (NT5_NaMN(1)).ET = 10
%      (NT5_NaMN(1)).scaling = 0.001
%      (NT5_NaMN(1)).kcat = 2.8
%      (NT5_NaMN(1)).Km = 3.5
%      (SIRT(1)).ET = 10
%      (SIRT(1)).scaling = 0.001
%      (SIRT(1)).Kcat = 0.67
%      (SIRT(1)).Km = 0.029
%      (SIRT(1)).Ki = 0.06
%      (NAD_consumption_without_Nam_inhibition(1)).ET = 10
%      (NAD_consumption_without_Nam_inhibition(1)).scaling = 0.001
%      (NAD_consumption_without_Nam_inhibition(1)).kcat = 1
%      (NAD_consumption_without_Nam_inhibition(1)).Km = 1
%      (NADA(1)).ET = 400
%      (NADA(1)).scaling = 0.001
%      (NADA(1)).Kcat = 0.69
%      (NADA(1)).Km = 0.009
%      (NADA(1)).Ki = 0.1
%      (NAMPT(1)).ET = 0
%      (NAMPT(1)).scaling = 0.001
%      (NAMPT(1)).Ki = 0.0021
%      (NAD_binding(1)).k1 = 100
%      (NAD_binding(1)).k2 = 10
%      (NRK1_NaMN(1)).ET = 10
%      (NRK1_NaMN(1)).scaling = 0.001
%      (NRK1_NaMN(1)).kcat = 0.23
%      (NRK1_NaMN(1)).Km = 0.0034
%      (PNP_NAR(1)).E_T = 10
%      (PNP_NAR(1)).turnover = 40
%      (PNP_NAR(1)).scaling = 0.001
%      (PNP_NAR(1)).Km = 1.48
%      (NRK1_NMN(1)).ET = 10
%      (NRK1_NMN(1)).scaling = 0.001
%      (NRK1_NMN(1)).kcat = 0.23
%      (NRK1_NMN(1)).Km = 0.0034
     v = 1e-05;
     Cytosol = 1;
     NamPT_compartment = 0.001;
%      NADA_compartment = 0.001
     celldivision_rate_NamPT_compartment = 2.8e-05;
%      celldivision_rate_NADA_compartment = 2.8e_05
     

% MV = NamPT_compartment .* ((SIRT_ET .* (SIRT_scaling .* (SIRT_Kcat .* x(:,9) .* H3_ac ./ (((SIRT_Km + x(:,9)) .* (1 + x(:,8)./(SIRT_Ki))))))))


% 
% %%      Initial Conditions:
%      x0(1)=1 ;           
% %      Nam_ex=1
%      x0(2)=1 ;          
% %      NADbound = 1;
%      x0(3)=1   ;        
% %      NR = 1;;
%      x0(4)=0   ;         
% %      NaAD = 0;
%      x0(5)=0  ;          
% %      NMN = 0;
%      x0(6)=0  ;          
% %      NA = 0;
%      x0(7)=0   ;         
% %      NaMN = 0;
%      x0(8)=0  ;          
% %      NAM = 0;
%      x0(9)= 0  ;         
% %      NAD = 0;
%      x0(10)= 0 ;        
% %      NAR = 0;


%      NADA_compartment.NADbound = 1;
%      NADA_compartment.NR = 1;
%      NADA_compartment.NaAD = 0;
%      NADA_compartment.NMN = 0;
%      NADA_compartment.NA = 0;
%      NADA_compartment.NaMN = 0;
%      NADA_compartment.NAM = 0;
%      NADA_compartment.NAD = 0;
%      NADA_compartment.NAR = 0;;
     Gln = 1;
     DNA_damage = 1;
     DNA_ADPR = 1;
     H3_ac = 1;
     H3_deac = 1;
     Glu = 1;
     Pi = 1;
     methyl_NAM = 0;
     H2O = 1;
     PPi = 0.012;
     AMP = 1;
     ADP = 1;
     SAM = 0.08;
     SAH = 1;
     PRPP = 1;
     ATP = 1;




      NAD_efflux = NamPT_compartment * celldivision_rate_NamPT_compartment * NAD;
%     (NAD_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NAD
      Nam_uptake = 0.001*v;
     
%      Repeated Assignments:
%     celldivision_rate_NADA_compartment = (Initial for cell devision all)
      celldivision_rate_NamPT_compartment = (Initial_for_cell_devision_all);


     NMN_efflux = NamPT_compartment * celldivision_rate_NamPT_compartment * NMN;
     NAPRT = NamPT_compartment * (NAPRT_E_T * NAPRT_turnover * NA * ATP * NAPRT_scaling / NAPRT_Km + NA);
     NMNAT1_NMN = NamPT_compartment * (ETNMNAT * NMNAT1_NMN_scaling * ((NMNAT1_NMN_kcat_A * NMN / (NMNAT1_NMN_Km_A) - (NMNAT1_NMN_kcat_PA * NAD / (NMNAT1_NMN_Km_PA)) / (1+NMN / (NMNAT1_NMN_Km_A + NaMN / (NMNAT1_NMN_Km_B + NAD / (NMNAT1_NMN_Km_PA + NaAD / (NMNAT1_NMN_Km_PB))))))));
     NMNAT1_NaMN = NamPT_compartment * (ETNMNAT * NMNAT1_NaMN_scaling * ((NMNAT1_NaMN_kcat_A * NaMN / (NMNAT1_NaMN_Km_A - (NMNAT1_NaMN_kcat_PA * NaAD / (NMNAT1_NaMN_Km_PA) / (1+NaMN / (NMNAT1_NaMN_Km_A + NMN / (NMNAT1_NaMN_Km_B + NaAD / (NMNAT1_NaMN_Km_PA + NAD / (NMNAT1_NaMN_Km_PB))))))))));
     NADS = NamPT_compartment * ((NADS_E_T * (NADS_turnover * NaAD * ATP * (NADS_scaling / ((NADS_Km + NaAD))))));
     NT5_NMN = NamPT_compartment * ((NT5_NMN_ET * (NT5_NMN_scaling * (NT5_NMN_kcat * NMN / ((NT5_NMN_Km + NMN))))));
     PNP_NR = NamPT_compartment * ((PNP_NR_E_T * (PNP_NR_turnover * NR * Pi * (PNP_NR_scaling / ((PNP_NR_Km + NR))))));
     NNMT = NamPT_compartment * ((NNMT_ET * (NNMT_scaling * (NNMT_Kcat * SAM * NAM / (Kma * (Kmb + NAM) * (1 + (methyl_NAM)/(NNMT_Ki) + SAM*Kmb + SAM*NAM))))));
     Nam_efflux = NamPT_compartment * celldivision_rate_NamPT_compartment * NAM;
     NT5_NaMN = NamPT_compartment * ((NT5_NaMN_ET * (NT5_NaMN_scaling * (NT5_NaMN_kcat * NaMN / ((NT5_NaMN_Km + NaMN))))));
     SIRT = NamPT_compartment *((SIRT_ET * (SIRT_scaling * (SIRT_Kcat * (0.01/1.3)*NADT * H3_ac/(((SIRT_Km + (0.01/1.3)*NADT) * (1+NAM/(SIRT_Ki))))))));
     NAD_consumption_without_Nam_inhibition = NamPT_compartment*((NAD_consumption_without_Nam_inhibition_ET*(NAD_consumption_without_Nam_inhibition_scaling*(NAD_consumption_without_Nam_inhibition_kcat*NAD/((NAD_consumption_without_Nam_inhibition_Km+NAD))))));
     NADA = NamPT_compartment*((NADA_ET*(NADA_scaling*(NADA_Kcat*NAM/((NADA_Km+NAM+(NADA_Km*NA/(NADA_Ki))))))));
     NAMimport = NamPT_compartment*((Nam_transporter)*Nam_ex-(Nam_transporter)*NAM);
     % pAMPK will increase NAMPT
     NAMPT = NamPT_compartment*((NAMPT_ET*(NAMPT_scaling*(Kcat_NamPRT)*NAM*(0.8+0.2*pAMPK/15)/((Km_NamPRT)+NAM+(Km_NamPRT)*NAD/(NAMPT_Ki)))));
     NAR_efflux = NamPT_compartment*celldivision_rate_NamPT_compartment*NAR;
     NA_efflux = NamPT_compartment*celldivision_rate_NamPT_compartment*NA;
     NAD_binding = NamPT_compartment*((NAD_binding_k1*NAD-(NAD_binding_k2*NADbound)));
     NRK1_NaMN = NamPT_compartment*((NRK1_NaMN_ET*(NRK1_NaMN_scaling*(NRK1_NaMN_kcat*NAR/((NRK1_NaMN_Km+NAR))))));
     PNP_NAR = NamPT_compartment*((PNP_NAR_E_T*(PNP_NAR_turnover*NAR*Pi*(PNP_NAR_scaling/((PNP_NAR_Km+NAR))))));
     NR_efflux = NamPT_compartment*celldivision_rate_NamPT_compartment*NR;
     NRK1_NMN = NamPT_compartment*((NRK1_NMN_ET*(NRK1_NMN_scaling*(NRK1_NMN_kcat*NR/((NRK1_NMN_Km+NR))))));

   











%% d/dt
     dNam_ex = (1/Cytosol) * (-(NAMimport) + 0.9625*(NAMimport) + (Nam_uptake));
     dNADbound = 1/NamPT_compartment * ((NAD_binding));
     dNR= 1/NamPT_compartment * ((NT5_NMN) - (PNP_NR) - (NR_efflux) - (NRK1_NMN));
     dNaAD = 1/NamPT_compartment * ((NMNAT1_NaMN) - (NADS));
     dNMN = 1/NamPT_compartment * (-(NMN_efflux) - (NMNAT1_NMN) - (NT5_NMN) + (NAMPT) + (NRK1_NMN));
     dNA = 1/NamPT_compartment * (-(NAPRT) + (NADA) - (NA_efflux) + (PNP_NAR));
     dNaMN = 1/NamPT_compartment * ((NAPRT) - (NMNAT1_NaMN) - (NT5_NaMN) + (NRK1_NaMN));
     dNAM = 1/NamPT_compartment * ((PNP_NR) - (NNMT) - (Nam_efflux) + (SIRT) + (NAD_consumption_without_Nam_inhibition) - (NADA) + (NAMimport) - (NAMPT));
     dNAD = 1/NamPT_compartment * ((NMNAT1_NMN) + (NADS) - (SIRT) - (NAD_consumption_without_Nam_inhibition) - (NAD_binding) - (NAD_efflux));
     dNAR = 1/NamPT_compartment * ((NT5_NaMN) - (NAR_efflux) - (NRK1_NaMN) - (PNP_NAR));
     dSIRT = 1000*(1000*SIRT-mySIRT);
%      d(NADA_compartment.NADbound)/dt = 1/NADA_compartment*((NAD_binding(1)))
%      d(NADA_compartment.NR)/dt = 1/NADA_compartment*((NT5_NMN(1)) _ (PNP_NR(1)) _ (NR_efflux(1)) _ (NRK1_NMN(1)))
%      d(NADA_compartment.NaAD)/dt = 1/NADA_compartment*((NMNAT1_NaMN(1)) _ (NADS(1)))
%      d(NADA_compartment.NMN)/dt = 1/NADA_compartment*(_(NMN_efflux(1)) _ (NMNAT1_NMN(1)) _ (NT5_NMN(1)) + (NAMPT(1)) + (NRK1_NMN(1)))
%      d(NADA_compartment.NA)/dt = 1/NADA_compartment*(_(NAPRT(1)) + (NADA(1)) _ (NA_efflux(1)) + (PNP_NAR(1)))
%      d(NADA_compartment.NaMN)/dt = 1/NADA_compartment*((NAPRT(1)) _ (NMNAT1_NaMN(1)) _ (NT5_NaMN(1)) + (NRK1_NaMN(1)))
%      d(NADA_compartment.NAM)/dt = 1/NADA_compartment*((PNP_NR(1)) _ (Nam_efflux(1)) + (SIRT(1)) + (NAD_consumption without Nam inhibition(1)) _ (NADA(1)) + (NAM import(1)) _ (NAMPT(1)))
%      d(NADA_compartment.NAD)/dt = 1/NADA_compartment*((NMNAT1_NMN(1)) + (NADS(1)) _ (SIRT(1)) _ (NAD_consumption without Nam inhibition(1)) _ (NAD_binding(1)) _ (NAD_efflux(1)))
%      d(NADA_compartment.NAR)/dt = 1/NADA_compartment*((NT5_NaMN(1)) _ (NAR_efflux(1)) _ (NRK1_NaMN(1)) _ (PNP_NAR(1)))
     
%  NamPT_compartment * ((SIRT_ET * (SIRT_scaling * (SIRT_Kcat * NAD * H3_ac/(((SIRT_Km + NAD) * (1+NAM/(SIRT_Ki))))))))

T=[dNam_ex; dNADbound; dNR; dNaAD; dNMN; dNA; dNaMN; dNAM; dNAD; dNAR; dSIRT]*60;
%  Convert rates from 1/s to 1/min



% %Fluxes:


%      (NMN_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NMN
%      (NAPRT(1)) = NADA_compartment*((NAPRT(1)).E_T*(NAPRT(1)).turnover*NADA_compartment.NA*NADA_compartment.ATP*(NAPRT(1)).scaling/((NAPRT(1)).Km+NADA_compartment.NA))
%      (NMNAT1_NMN(1)) = NADA_compartment*((E_T NMNAT)*(NMNAT1_NMN(1)).scaling*((NMNAT1_NMN(1)).kcat_A*NADA_compartment.NMN/(NMNAT1_NMN(1)).Km_A_(NMNAT1_NMN(1)).kcat_PA*NADA_compartment.NAD/(NMNAT1_NMN(1)).Km_PA)/(1+NADA_compartment.NMN/(NMNAT1_NMN(1)).Km_A+NADA_compartment.NaMN/(NMNAT1_NMN(1)).Km_B+NADA_compartment.NAD/(NMNAT1_NMN(1))._PA+NADA_compartment.NaAD/(NMNAT1_NMN(1)).Km_PB))
%      (NMNAT1_NaMN(1)) = NADA_compartment*((E_T NMNAT)*(NMNAT1_NaMN(1)).scaling*((NMNAT1_NaMN(1)).kcat_A*NADA_compartment.NaMN/(NMNAT1_NaMN(1)).Km_A_(NMNAT1_NaMN(1)).kcat_PA*NADA_compartment.NaAD/(NMNAT1_NaMN(1)).Km_PA)/(1+NADA_compartment.NaMN/(NMNAT1_NaMN(1)).Km_A+NADA_compartment.NMN/(NMNAT1_NaMN(1)).Km_B+NADA_compartment.NaAD/(NMNAT1_NaMN(1)).Km_PA+NADA_compartment.NAD/(NMNAT1_NaMN(1)).Km_PB))
%      (NADS(1)) = NADA_compartment*((NADS(1)).E_T*(NADS(1)).turnover*NADA_compartment.NaAD*NADA_compartment.ATP*(NADS(1)).scaling/((NADS(1)).Km+NADA_compartment.NaAD))
%      (NT5_NMN(1)) = NADA_compartment*((NT5_NMN(1)).ET*(NT5_NMN(1)).scaling*(NT5_NMN(1)).kcat*NADA_compartment.NMN/((NT5_NMN(1)).Km+NADA_compartment.NMN))
%      (PNP_NR(1)) = NADA_compartment*((PNP_NR(1)).E_T*(PNP_NR(1)).turnover*NADA_compartment.NR*NADA_compartment.Pi*(PNP_NR(1)).scaling/((PNP_NR(1)).Km+NADA_compartment.NR))
%      (Nam_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NAM
%      (NT5_NaMN(1)) = NADA_compartment*((NT5_NaMN(1)).ET*(NT5_NaMN(1)).scaling*(NT5_NaMN(1)).kcat*NADA_compartment.NaMN/((NT5_NaMN(1)).Km+NADA_compartment.NaMN))
%      (SIRT(1)) = NADA_compartment*((SIRT(1)).ET*(SIRT(1)).scaling*(SIRT(1)).Kcat*NADA_compartment.NAD*NADA_compartment.H3_ac/(((SIRT(1)).Km+NADA_compartment.NAD)*(1+NADA_compartment.NAM/(SIRT(1)).Ki)))
%      (NAD_consumption without Nam inhibition(1)) = NADA_compartment*((NAD_consumption without Nam inhibition(1)).ET*(NAD_consumption without Nam inhibition(1)).scaling*(NAD_consumption without Nam inhibition(1)).kcat*NADA_compartment.NAD/((NAD_consumption without Nam inhibition(1)).Km+NADA_compartment.NAD))
%      (NADA(1)) = NADA_compartment*((NADA(1)).ET*(NADA(1)).scaling*(NADA(1)).Kcat*NADA_compartment.NAM/((NADA(1)).Km+NADA_compartment.NAM+(NADA(1)).Km*NADA_compartment.NA/(NADA(1)).Ki))
%      (NAM import(1)) = NADA_compartment*((Nam transporter)*Nam_ex_(Nam transporter)*NADA_compartment.NAM)
%      (NAMPT(1)) = NADA_compartment*((NAMPT(1)).ET*(NAMPT(1)).scaling*(Kcat NamPRT)*NADA_compartment.NAM/((Km NamPRT)+NADA_compartment.NAM+(Km NamPRT)*NADA_compartment.NAD/(NAMPT(1)).Ki))
%      (NAR_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NAR
%      (NA_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NA
%      (NAD_binding(1)) = NADA_compartment*((NAD_binding(1)).k1*NADA_compartment.NAD_(NAD_binding(1)).k2*NADA_compartment.NADbound)
%      (NRK1_NaMN(1)) = NADA_compartment*((NRK1_NaMN(1)).ET*(NRK1_NaMN(1)).scaling*(NRK1_NaMN(1)).kcat*NADA_compartment.NAR/((NRK1_NaMN(1)).Km+NADA_compartment.NAR))
%      (PNP_NAR(1)) = NADA_compartment*((PNP_NAR(1)).E_T*(PNP_NAR(1)).turnover*NADA_compartment.NAR*NADA_compartment.Pi*(PNP_NAR(1)).scaling/((PNP_NAR(1)).Km+NADA_compartment.NAR))
%      (NR_efflux(1)) = NADA_compartment*celldivision_rate_NADA_compartment*NADA_compartment.NR
%      (NRK1_NMN(1)) = NADA_compartment*((NRK1_NMN(1)).ET*(NRK1_NMN(1)).scaling*(NRK1_NMN(1)).kcat*NADA_compartment.NR/((NRK1_NMN(1)).Km+NADA_compartment.NR))

%      NADA_compartment.Gln = 1;
%      NADA_compartment.DNA_damage = 1;
%      NADA_compartment.DNA_ADPR = 1
%      NADA_compartment.H3_ac = 1;
%      NADA_compartment.H3_deac = 1;
%      NADA_compartment.Glu = 1;
%      NADA_compartment.Pi = 1;
%      NADA_compartment.H2O = 1;
%      NADA_compartment.PPi = 0.012;
%      NADA_compartment.AMP = 1;
%      NADA_compartment.ADP = 1;
%      NADA_compartment.PRPP = 1;
%      NADA_compartment.ATP = 1;




end


function yprime=circadian_mtorc(t,y,mtorc1,ampk,PP,aging,stac_params,nad_params)
%% Notations
% y(1)=per;
% y(2)=cry;
% y(3)=rev;
% y(4)=ror;
% y(5)=bmal;
% y(6)=Prot_per;
% y(7)=Prot_cry;
% y(8)=Prot_rev;
% y(9)=Prot_ror; 
% y(10)=Prot_bmal;
% y(11)=PC;
% y(12)=CB;
% y(13)=nampt;
% y(14)=Prot_nampt;
% y(15)=NAD;
% y(16)=dbp; 

% simulate STAC and NAD
% params = [amplitude, timing, duration, rise time]
stac_scaling = 1-stac_params(1)*get_P(t,stac_params(2),stac_params(3),stac_params(4),24);
nad_scaling  = 1+nad_params(1)*get_P(t,nad_params(2),nad_params(3),nad_params(4),24);

% k_CBM=0.0242;

per=y(1);
cry=y(2);
rev=y(3);
ror=y(4);
bmal=y(5);
Prot_per=y(6);
Prot_cry=y(7);
Prot_rev=y(8);
Prot_ror=y(9);
Prot_bmal=y(10);
PC=y(11);
CB=y(12);
nampt=y(13);
Prot_nampt=y(14);
NAD=y(15);
dbp=y(16);

K_bm=PP(1);
kp_bmal=PP(2);
dp_bmal=PP(3);
kass_cb=PP(4);
kdiss_cb=PP(5);
k_CBM=PP(6);
d_cb=PP(7); 
K_BinhibitM=PP(8);
K_PM=PP(9);

kass_pc=PP(10);
kdiss_pc=PP(11);
dp_cry=PP(12);
dp_per=PP(13);
kp_cry=PP(14);
kp_per=PP(15); 
kp_rev=PP(16);
kp_ror=PP(17);
dm_per=PP(18);
Vmax_per=PP(19);


%% Define constant

% Table S3A (unit: h^-1)  mRNA and protein degradatation rates

% MANUAL FIT
dm_bmal=0.827333442085 * 1.05;    % Bmal mRNA degradation rate

% MANUAL FIT
dm_cry=0.319848706181 * 0.9;     % Cry mRNA degradation rate

dm_dbp=0.384486532062;     % Dbp mRNA degradation rate

dm_nampt=0.814311309051;   % Nampt mRNA degradation rate

% dm_per=0.323114598647;     % Per mRNA degradation rate

dm_rev=4.0479072846;       % Rev-Erb mRNA degradation rate

dm_ror=0.246575760727;     % Ror mRNA degradation rate

% dp_bmal=0.186413466797;    % BMAL protein degradation rate

% dp_cry=0.599026119971;     % CRY protein degradation rate

dp_nampt=49.8841023982;    % NAMPT protein degradation rate

% dp_per=10.9446515392;      % PER protein degradation rate

dp_rev=0.281564063057;     % REV-ERB protein degradation rate

dp_ror=0.0340112196281 * 1.6;    %  ROR degradation rate

% d_cb=0.197714012552;       % CLOCK-BMAL complex degradation rate
% d_cb=0.0988;

d_pc=0.609290138329;       % PER-CRY complex degradation rate








% Table S3B Complexation kinetic rates


% kass_cb=0.0162923358655;    % (nmol^-1.l.h^-1) CLOCK-BMAL association rate

% kass_pc=12.302005485;       % (nmol^-1.l.h^-1) PER-CRY association rate

% kdiss_cb=0.00613502224231;  % (h^-1) CLOCK-BMAL dissociation rate

% kdiss_pc=0.0365867175408;       % (h^-1) PER-CRY dissociation rate









% Table S3C  (unit: nmol.l^-1.h^-1)   Maximal transcription rates


Vmax_bmal=0.0916862770193;  % Bmal maximal transcription rate

Vmax_cry=0.702216664807;    % Cry maximal transcription rate

Vmax_dbp=0.0802009609453;   % Dbp maximal transcription rate

Vmax_nampt=3.49035201578;   % Nampt maximal transcription rate

% Vmax_per= 0.730201742662;   % Per maximal transcription rate

Vmax_rev= 1.12297601784;    % Rev-Erb maximal transcription rate

Vmax_ror= 6.9843472736;     % Ror maximal transcription rate








% Table S3D (dimensionless) Activation ratios


fold_bmal=15.9258093373;    % activation ratio of Bmal by ROR

fold_cry= 1.1604489571;     % activation ratio of Cry by CLOCK-BMAL

fold_dbp=400.0;             % activation ratio of Dbp by CLOCK-BMAL

fold_nampt=1.57880681573;   % activation ratio of Nampt by CLOCK-BMAL

fold_per=12.977351123;      % activation ratio of Per by CLOCK-BMAL

fold_rev=73.2342431701;     % activation ratio of Rev-Erb by CLOCK-BMAL

fold_ror=335.923333883;     % activation ratio of Ror by CLOCK-BMAL







% Table S3E (nmol^-1.l) Regulation threshold

Ka_bmal_ror=0.00560038941378; % Regulation threshold of Bmal by ROR

Ka_cry_cb=1.0089387144;       % Regulation threshold of Cry by CLOCK-BMAL

Ka_dbp_cb=0.308745016237;     % Regulation threshold of Dbp by CLOCK-BMAL

Ka_nampt_cb=3.54586790835;    % Regulation threshold of Nampt by CLOCK-BMAL

Ka_per_cb=2.03485134763;      % Regulation threshold of Per by CLOCK-BMAL

Ka_rev_cb=0.260846828116;     % Regulation threshold of Rev-Erb by CLOCK-BMAL

Ka_ror_cb=0.266407416327;     % Regulation threshold of Ror by CLOCK-BMAL

Ki_bmal_rev0=0.0108449480001; % Regulation threshold of Bmal by REV-ERB

Ki_cry_rev0=0.248955507809;   % Regulation threshold of Cry by REV-ERB

Crev=0;

Ki_cry_pc=0.00338463577329;   % Regulation threshold of Cry by PER-CRY

Ki_dbp_pc=2.23913672671;      % Regulation threshold of Dbp by PER-CRY

Ki_nampt_pc=0.0137106537972;  % Regulation threshold of Nampt by PER-CRY

Ki_per_pc=0.273493946059;     % Regulation threshold of Per by PER-CRY

Ki_rev_pc=28.5885406354;      % Regulation threshold of Rev-Erb by PER-CRY

Ki_ror_pc=0.0072858432208;    % Regulation threshold of Ror by PER-CRY








% Table S3F (dimensionless) Hill coefficient

hill_bmal_rev=4.32985205032;   % Hill coeff., regulation of Bmal by REV-ERB

hill_bmal_ror=1.83992599578;   % Hill coeff., regulation of Bmal by ROR

hill_cry_cb=9.1109447538;      % Hill coeff., regulation of Cry by CLOCK-BMAL

hill_cry_pc=2.43715119318;     % Hill coeff., regulation of Cry by PER-CRY

hill_cry_rev=4.20952050286;    % Hill coeff., regulation of Cry by REV-ERB

hill_dbp_cb=7.32066818222;     % Hill coeff., regulation of Dbp by CLOCK-BMAL

hill_dbp_pc=10.4312927466;     % Hill coeff., regulation of Dbp by PER-CRY

hill_nampt_cb=1.91470474775;   % Hill coeff., regulation of Nampt by CLOCK-BMAL

hill_nampt_pc=1.34080593157;   % Hill coeff., regulation of Nampt by PER-CRY

hill_per_cb=8.52414053707;     % Hill coeff., regulation of Per by CLOCK-BMAL

hill_per_pc=8.53897990872;     % Hill coeff., regulation of Per by PER-CRY

hill_rev_cb=9.83701536835;     % Hill coeff., regulation of Rev-Erb by CLOCKBMAL

hill_rev_pc=3.31257899336;     % Hill coeff., regulation of Rev-Erb by PER-CRY

hill_ror_cb=9.36456505724;     % Hill coeff., regulation of Ror by CLOCK-BMAL

hill_ror_pc=1.84102439743;     % Hill coeff., regulation of Ror by PER-CRY









% Table S3G (unit: molecules per hour per mRNA) Translation rates


% kp_bmal=0.628507384997;      % Bmal translation rate

% kp_cry=3.76912711677;        % Cry translation rate

kp_nampt=58.9647983852;      % Nampt translation rate

% kp_per=13.2909782781;        % Per translation rate

% kp_rev=0.0513221194724;      % Rev-Erb translation rate

% kp_ror=0.0412765888526;      % Ror translation rate







% Table S3H (dimensionless) Protein stability modulation constants

m_cry_ampk=0.07940386211;   % modulation of CRY stability by AMPK

m_nampt_ampk=0.6352243791;  % modulation of NAMPT stability by AMPK

m_per_ampk=0.005243953731;  % modulation of PER stability by AMPK

m_per_sirt=0.005452322816;  % modulation of PER stability by SIRT







% Table S3I  NAD kinetics, Sirt1 and PGC1a activity

Vsirt=0.915854846427 * aging(2);      % (dimensionless) Maximum SIRT1 activity

Ksirt=0.75;                % (nmol.l^-1)  Value of [NAD] at which SIRT1 activity is half of maximum

d_nad=378.009130749;       % (h^-1) Rate of transformation of NAD into NAM

Knad=0.321746039086;       % (nmol.l^-1) Value of [NAD] at which transformation into NAM is at half of maximum rate

NAD_basal=0.9116166306 * aging(1);    % (nmol.l^-1) Value of [NAD] below which transformation of NAD into NAM is inactivated

Vnad=296.3933869;          % (molecule per hour per NAMPT protein) Maximum regeneration rate of NAD

NAD_tot=4.166901679 * aging(1);       % (nmol.l^-1) Total concentration of NAD and NAM
NAD_tot=NAD_tot * nad_scaling;

Knam =2.76496;             % (nmol.l^-1) Value of NAM at which NAD salvage rate is half of maximum

Vpg=24.06372088;           % (nmol^-1.h) Maximum activity of PGC1a

Kpg1=0.046630145542;       % (dimensionless) Michaelis-Menten constant for phosphorylation of PGC1a by AMPK

Kpg2=12.3526351747;        % (dimensionless) Michaelis-Menten constant for deacetylation of PGC1a by SIRT1





% Table S3J: Pulse parameters

tc1=4.38900149;            % (h) Timing of the first AMPK pulse
tc2=15.75;                 % (h) Timing of the second AMPK pulse
tc3=18.875;                % (h) Time of maximal nuclear PGC1a abundance
Td1=2.25;                  % (h) Duration of the first AMPK pulse
Td2=1.5;                   % (h) Duration of the first AMPK pulse
Td3=15.25;                 % (h) Duration of the nuclear PGC1a presence
Tr1=2.6;                   % (h) Rise time of the first AMPK pulse
Tr2=1.8;                   % (h) Rise time of the second AMPK pulse
Tr3=0.5;                   % (h) Rise time of nuclear PGC1a
amp1=6.0;                  % (N/A) Amplitude of the first AMPK pulse
amp2=0.9778008;            % (N/A) Amplitude of the second AMPK pulse
amp3=0.803062;             % (nmol^-1.l) Amplitude of the nuclear PGC1a abundance pulse





% Table S3K: Chronotherapy timings

tc4=13.664;                % Timing of the agonist pulse

Td4=2.83718;               % Duration of the agonist pulse

Tr4=1.86794;               % Rise time of the agonist pulse

amp4=0.465852;             % Amplitude of the agonist pulse





% Table S3L: Miscellaneous constants used to describe perturbations

Csirt=1;                  % 1 for WT, LKB1KO, HFD, fasting; 0 for SIRT1KO

Campk=1;                  % 1 for WT, SIRT1KO; 0.0375 for LKB1KO; 0.05 for HFD and fasting       

Cpgc1=1;                  % 1 for WT 

cpg=1.009872; 

offs=0.02;                % 0.02 for WT,  SIRT1KO, LKB1KO, HFD; 2.6 for fasting 



%% Derivatives 

% dper

Ksirt = Ksirt * stac_scaling;
Act_SIRT=(Csirt*Vsirt*NAD)/(Ksirt+NAD);
% Act_SIRT=1.5;
top_per=Vmax_per*(1+fold_per*(CB/(Ka_per_cb*(1+Act_SIRT)))^hill_per_cb);

bottom_per=1+(CB/(Ka_per_cb*(1+Act_SIRT)))^hill_per_cb*(1+(PC/Ki_per_pc)^hill_per_pc);

dper=-dm_per*per+top_per/bottom_per;



% dcry

top_cry=Vmax_cry*(1+fold_cry*(CB/(Ka_cry_cb*(1+Act_SIRT)))^hill_cry_cb);

bottom_cry1=1+(CB/(Ka_cry_cb*(1+Act_SIRT)))^hill_cry_cb*(1+(PC/Ki_cry_pc)^hill_cry_pc);

P4=get_P(t,tc4,Td4,Tr4,24);

agonist_rev=amp4*P4;

Ki_cry_rev=Ki_cry_rev0/(1+Crev*agonist_rev);


bottom_cry2=1+(Prot_rev/Ki_cry_rev)^hill_cry_rev; 

dcry=-dm_cry*cry+top_cry/bottom_cry1/bottom_cry2;





% drev

top_rev=Vmax_rev*(1+fold_rev*(CB/(Ka_rev_cb*(1+Act_SIRT)))^hill_rev_cb);

bottom_rev=1+(CB/(Ka_rev_cb*(1+Act_SIRT)))^hill_rev_cb*(1+(PC/Ki_rev_pc)^hill_rev_pc);

drev=-dm_rev*rev+top_rev/bottom_rev;




% dror

top_ror=Vmax_ror*(1+fold_ror*(CB/(Ka_ror_cb*(1+Act_SIRT)))^hill_ror_cb);

bottom_ror=1+(CB/(Ka_ror_cb*(1+Act_SIRT)))^hill_ror_cb*(1+(PC/Ki_ror_pc)^hill_ror_pc);

dror=-dm_ror*ror+top_ror/bottom_ror;




% dbmal

P1=get_P(t,tc1,Td1,Tr1,24);

P2=get_P(t,tc2,Td2,Tr2,24);

P3=get_P(t,tc3,Td3,Tr3,24);

Act_AMPK=Campk*(amp1*P1+amp2*P2)+(1-Campk)*offs;
% Act_AMPK=ampk;

Prot_PGC1a=amp3*P3;

Act_PGC1a=(Cpgc1*cpg*Vpg*Act_AMPK*Act_SIRT*Prot_PGC1a)/(1+Act_AMPK/Kpg1*(1+Act_SIRT/Kpg2));

Ki_bmal_rev=Ki_bmal_rev0/(1+Crev*agonist_rev);

top_bmal=Vmax_bmal*(1+fold_bmal*(1+Act_PGC1a)*(Prot_ror/Ka_bmal_ror)^hill_bmal_ror);

bottom_bmal=1+(Prot_rev/Ki_bmal_rev)^hill_bmal_rev+(Prot_ror/Ka_bmal_ror)^hill_bmal_ror;

dbmal=-dm_bmal*bmal+top_bmal/bottom_bmal;





% dProt_per
dProt_per=-dp_per*(1+m_per_sirt*Act_SIRT+m_per_ampk*Act_AMPK)*Prot_per...
    +kp_per*per...
    -(kass_pc*Prot_cry*Prot_per-kdiss_pc*PC); 




% dProt_cry
dProt_cry=-dp_cry*(1+m_cry_ampk*Act_AMPK)*Prot_cry...
    +kp_cry*cry...
    -(kass_pc*Prot_cry*Prot_per-kdiss_pc*PC);




% dProt_rev
dProt_rev=-dp_rev*Prot_rev+kp_rev*rev;



% dProt_ror
dProt_ror=-dp_ror*Prot_ror+kp_ror*ror;





% dProt_bmal
dProt_bmal=-dp_bmal*Prot_bmal...
    +kp_bmal*bmal...
    -(kass_cb*Prot_bmal-kdiss_cb*CB)...
    +bmal*mtorc1*K_bm;


% [dProt_bmal -dp_bmal*Prot_bmal,     +kp_bmal*bmal,     -(kass_cb*Prot_bmal-kdiss_cb*CB),     +bmal*mtorc1*K_bm]
% pause

%[bmal,mtorc1,K_bm],pause
%[kdiss_cb*CB, -kass_cb*Prot_bmal, bmal*mtorc1*K_bm, kp_bmal*bmal, -dp_bmal*Prot_bmal],pause


% dPC
dPC=(kass_pc*Prot_cry*Prot_per-kdiss_pc*PC)-d_pc*PC;



% dCB
dCB=(kass_cb*Prot_bmal-kdiss_cb*CB)-d_cb*CB - mtorc1*k_CBM*CB;

% [dCB, kass_cb*Prot_bmal, -kdiss_cb*CB, -d_cb, CB, -mtorc1, k_CBM*CB],
% pause



% dnampt
top_nampt=Vmax_nampt*(1+fold_nampt*(CB/(Ka_nampt_cb*(1+Act_SIRT)))^hill_nampt_cb);

bottom_nampt=1+(CB/(Ka_nampt_cb*(1+Act_SIRT)))^hill_nampt_cb*(1+(PC/Ki_nampt_pc)^hill_nampt_pc);

dnampt=-dm_nampt*nampt+top_nampt/bottom_nampt;





% dProt_nampt

top_Prot_nampt=dp_nampt*Prot_nampt;

bottom_Prot_nampt=1+m_nampt_ampk*Act_AMPK;

dProt_nampt=-top_Prot_nampt/bottom_Prot_nampt+kp_nampt*nampt;





% dNAD
top_NAD_1=d_nad*(NAD-NAD_basal);

bottom_NAD_1=Knad+NAD-NAD_basal;

top_NAD_2=Vnad*Prot_nampt*(NAD_tot-NAD);

bottom_NAD_2=Knam+NAD_tot-NAD; 

dNAD=-top_NAD_1/bottom_NAD_1+top_NAD_2/bottom_NAD_2;








% ddbp

top_dbp=Vmax_dbp*(1+fold_dbp*(CB/(Ka_dbp_cb*(1+Act_SIRT)))^hill_dbp_cb);

bottom_dbp=1+(CB/(Ka_dbp_cb*(1+Act_SIRT)))^hill_dbp_cb*(1+(PC/Ki_dbp_pc)^hill_dbp_pc);

ddbp=-dm_dbp*dbp+top_dbp/bottom_dbp;

 


yprime=[dper; dcry; drev; dror; dbmal; dProt_per; dProt_cry; dProt_rev; dProt_ror; dProt_bmal; dPC; dCB; dnampt; dProt_nampt; dNAD; ddbp];
% sirt_scaling = 0.0026/1.5;
sirt_scaling = 0.005/1.5;
yprime=[yprime; Act_SIRT*sirt_scaling];

% disp( ['time = ',num2str(t)])

end



function z=function_6(IRS_p,IRS_p_phos_by_p70_S6K_pT229_pT389,p70_S6K_pT229_pT389), z=(IRS_p_phos_by_p70_S6K_pT229_pT389*IRS_p*p70_S6K_pT229_pT389);end

function z=function_3(Amino_Acids,mTORC2,mTORC2_S2481_phos_by_Amino_Acids), z=(mTORC2_S2481_phos_by_Amino_Acids*mTORC2*Amino_Acids);end

function z=function_14(Akt_S473_phos_by_mTORC2_pS2481_second,Akt_pT308,mTORC2_pS2481), z=(Akt_S473_phos_by_mTORC2_pS2481_second*Akt_pT308*mTORC2_pS2481);end

function z=function_17(TSC1_TSC2,TSC1_TSC2_pS1387,mTORC1_pS2448,mTORC1_pS2448_dephos_by_TSC1_TSC2), z=(mTORC1_pS2448_dephos_by_TSC1_TSC2*mTORC1_pS2448*(TSC1_TSC2+TSC1_TSC2_pS1387));end

function z=function_20(PI3K_p_PDK1,p70_S6K,p70_S6K_T229_phos_by_PI3K_p_PDK1_first), z=(p70_S6K_T229_phos_by_PI3K_p_PDK1_first*p70_S6K*PI3K_p_PDK1);end

function z=function_19(PI3K_variant_p,mTORC2,mTORC2_S2481_phos_by_PI3K_variant_p), z=(mTORC2_S2481_phos_by_PI3K_variant_p*mTORC2*PI3K_variant_p);end

function z=function_2(AMPK,AMPK_T172_phos_by_Amino_Acids,Amino_Acids), z=(AMPK_T172_phos_by_Amino_Acids*AMPK*Amino_Acids);end

function z=function_4(IR_beta,IR_beta_phos_by_Insulin,Insulin), z=(IR_beta_phos_by_Insulin*IR_beta*Insulin);end

function z=function_23(mTORC1_pS2448,p70_S6K_T389_phos_by_mTORC1_pS2448_second,p70_S6K_pT229), z=(p70_S6K_T389_phos_by_mTORC1_pS2448_second*p70_S6K_pT229*mTORC1_pS2448);end

function z=function_22(PI3K_p_PDK1,p70_S6K_T229_phos_by_PI3K_p_PDK1_second,p70_S6K_pT389), z=(p70_S6K_T229_phos_by_PI3K_p_PDK1_second*p70_S6K_pT389*PI3K_p_PDK1);end

function z=function_16(Akt_pT308,Akt_pT308_pS473,TSC1_TSC2,TSC1_TSC2_T1462_phos_by_Akt_pT308), z=(TSC1_TSC2_T1462_phos_by_Akt_pT308*TSC1_TSC2*(Akt_pT308+Akt_pT308_pS473));end

function z=function_11(Akt,Akt_T308_phos_by_PI3K_p_PDK1_first,PI3K_p_PDK1), z=(Akt_T308_phos_by_PI3K_p_PDK1_first*Akt*PI3K_p_PDK1);end

function z=function_27(PRAS40_S183_phos_by_mTORC1_pS2448_second,PRAS40_pT246,mTORC1_pS2448), z=(PRAS40_S183_phos_by_mTORC1_pS2448_second*PRAS40_pT246*mTORC1_pS2448);end

function z=function_26(Akt_pT308,Akt_pT308_pS473,PRAS40_T246_phos_by_Akt_pT308_second,PRAS40_pS183), z=(PRAS40_T246_phos_by_Akt_pT308_second*PRAS40_pS183*(Akt_pT308+Akt_pT308_pS473));end

function z=function_7(IRS,IRS_phos_by_p70_S6K_pT229_pT389,p70_S6K_pT229_pT389), z=(IRS_phos_by_p70_S6K_pT229_pT389*IRS*p70_S6K_pT229_pT389);end

function z=function_15(AMPK_pT172,TSC1_TSC2,TSC1_TSC2_S1387_phos_by_AMPK_pT172), z=(TSC1_TSC2_S1387_phos_by_AMPK_pT172*TSC1_TSC2*AMPK_pT172);end

function z=function_9(IR_beta_pY1146,PI3K_variant,PI3K_variant_phos_by_IR_beta_pY1146), z=(PI3K_variant_phos_by_IR_beta_pY1146*PI3K_variant*IR_beta_pY1146);end

function z=function_12(Akt,Akt_S473_phos_by_mTORC2_pS2481_first,mTORC2_pS2481), z=(Akt_S473_phos_by_mTORC2_pS2481_first*Akt*mTORC2_pS2481);end

function z=function_18(Amino_Acids,mTORC1,mTORC1_S2448_activation_by_Amino_Acids), z=(mTORC1_S2448_activation_by_Amino_Acids*mTORC1*Amino_Acids);end

function z=function_10(AMPK,AMPK_T172_phos,IRS_p), z=(AMPK_T172_phos*AMPK*IRS_p);end

function z=function_13(Akt_T308_phos_by_PI3K_p_PDK1_second,Akt_pS473,PI3K_p_PDK1), z=(Akt_T308_phos_by_PI3K_p_PDK1_second*Akt_pS473*PI3K_p_PDK1);end

function z=function_8(IRS_p,PI3K_PDK1,PI3K_PDK1_phos_by_IRS_p), z=(PI3K_PDK1_phos_by_IRS_p*PI3K_PDK1*IRS_p);end

function z=function_21(mTORC1_pS2448,p70_S6K,p70_S6K_T389_phos_by_mTORC1_pS2448_first), z=(p70_S6K_T389_phos_by_mTORC1_pS2448_first*p70_S6K*mTORC1_pS2448);end

function z=function_5(IRS,IRS_phos_by_IR_beta_pY1146,IR_beta_pY1146), z=(IRS_phos_by_IR_beta_pY1146*IRS*IR_beta_pY1146);end

function z=function_24(PRAS40,PRAS40_S183_phos_by_mTORC1_pS2448_first,mTORC1_pS2448), z=(PRAS40_S183_phos_by_mTORC1_pS2448_first*PRAS40*mTORC1_pS2448);end

function z=function_25(Akt_pT308,Akt_pT308_pS473,PRAS40,PRAS40_T246_phos_by_Akt_pT308_first), z=(PRAS40_T246_phos_by_Akt_pT308_first*PRAS40*(Akt_pT308+Akt_pT308_pS473));end

function z=function_1(Amino_Acids,IRS,IRS_phos_by_Amino_Acids), z=(IRS_phos_by_Amino_Acids*IRS*Amino_Acids);end

% adding few functions representing operators used in SBML but not present directly 
% in either matlab or octave. 
% function z=pow(x,y),z=x^y;end
% function z=root(x,y),z=y^(1/x);end
% function z = piecewise(varargin)
% 	numArgs = nargin;
% 	result = 0;
% 	foundResult = 0;
% 	for k=1:2: numArgs-1
% 		if varargin{k+1} == 1
% 			result = varargin{k};
% 			foundResult = 1;
% 			break;
% 		end
% 	end
% 	if foundResult == 0
% 		result = varargin{numArgs};
% 	end
% 	z = result;
% end

