package gov.nih.ncats.ifx.qhts.utils;

public class HillConstants {

	
	//fitting
    public static double MIN_Y_RANGE = 20.0; //default
    public static double Y0 = 10.0;
    public static double PS_MIN = 20.0;  //to be removed
    public static double PS_MAX = 30.0; //not used
    public static double PI_MIN = -30.0;  //not used
    public static double PI_MAX = -20.0; //to be removed

    public static double SUPER_Y = 30.0;  //do not check range

    public static int PARTIAL_FIT_NO = 5;
    public static int PARTIAL_FIT_MASK_NO = 2;

    public static boolean CURVE_CLASS5 = false;
    public static double R2 = 0.30;
    public static double MIN_SLOPE = 0.1;
    public static double MAX_SLOPE = 5.0;

    public static double Y0_INF_COEF = 1.2; //limit of y0 and yinf

    public static String FIT_TYPE = "P4";

    public static double EC50_GRID_SIZE = 0.1;


    //fitting parameters
    public static String FIXED_Y0;
    public static String FIXED_YINF;
    public static String FIXED_SLOPE;
    public static String FIXED_AC50;
	public static String MIN_YINF;
	public static String MAX_YINF;
    public static boolean Y0_LESS_THAN_YINF = false;

    public static double MIN_LOG_EC1 = -10;

    public static boolean UPDATE_FIT_DATA = false;

    public static boolean YINF_FROM_DATA = false; //yinf from data
    public static boolean YINF_LESS_100 = false; //yinf <= 100%
    public static boolean YINF_MORE_80 = false; //yinf >  80%

	public static double OCCAM_YMAX = 100.0;
	
	//super active
	public static boolean AUTO_SUPER_ACTIVE = false;

    //biphase
    public static boolean WEIGHT_FOR_BIPHASE = true;

    //masking
    public static double P0 = 20.0; //% activity allowed for first data point

    public static double PMAX = -250.0; // max % activity; not used

    public static double THETA = 120.0;

    public static double V_DEPTH = 20.0;

    public static double TP0 = 20.0; //mask delta cutoff for first point
    public static double TPV = 30.0; //mask delta cutoff for v shape
    public static double TP = 70.0; //mask delta cutoff
    public static double TPN = 150.0; //mask delta cutoff for last point
    public static double TR = 5.0; //(y-mean)/std

    public static int MASK = 6; //frequency

    //Ruili uses "true"
    public static boolean BELL_MASK = true; //data point lower than peak
    //public static boolean BELL_MASK = false; //data point lower than peak
    
    public static boolean BELL_CHECK_REGRESSION = true; //data point lower than peak
    public static int BELL_CHECK_START = 2; // 1/2 of data points
	public static boolean BELL_ABS = true;

    public static boolean MASK_FLAG = true;

	public static String USER_MASK = "";

    //confidence interval
    public static boolean CI_FLAG = false;

    public static String USE_MASKED_POINT_OPTION = "No for Bell";

    //weighting
    public static String WEIGHT_SCHEME = "NONE";

    //classification
    public static double CLASSIFICATION_SD = 10.0;
    public static double CLASSIFICATION_SD_FACTOR = 3;
    public static double R2_CLASS1_CUTOFF = 0.90;
    public static boolean RELATIVE_ACTIVITY = true;
    public static boolean SHIFT_ACTIVITY = false;
    public static String CLASS1_CUTOFF = "80%";
    public static double RANGE_CUTOFF = 0.0;
	public static double SUPER_ACTIVE_RATIO = 0.8;
	public static double ROBUST_PCT = 80.0;
	public static double C124_FACTOR = 2.0;
	
}
