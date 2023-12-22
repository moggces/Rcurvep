package gov.nih.ncats.ifx.qhts.utils;

public class NCGCConstants {
	
	//fitting
    public static boolean REFIT = false;

    //multi-view correction
    public static double MIN_RESPONSE_TO_DRAW = 10.0;
    public static double COLOR_FACTOR = 2.0;

    //marvin
    public static int   MOL_PER_ROW = 3;

    //charting
    public static boolean CURVE_CLASS_RANGE = true;
    public static boolean AUTO_RANGE = true;
    public static double MIN_RANGE = -120.0;
    public static double MAX_RANGE = 120.0;
	public static boolean FIRST_CMPD_ONLY = true;
	public static boolean SHOW_CC4_FIT = false;

    //loading
    public static boolean ROLLBACK_ON_MISSING_PLATE_LOG = false;
    public static String INCLUDE_LAYER_INDEXES = null;

    //data
    public static boolean FILL_MISSING_DATA = true;
    public static boolean RAW_DATA = false; //if true, normalize to 100
    public static boolean YINF_FROM_DATA = false; //yinf from data
    public static boolean YINF_LESS_100 = true; //yinf <= 100%
    public static boolean Y0_AS_100 = false;



    //hill curve legend
    public static String LEGEND_SHOW_OPTION = "Sample ID";
    public static String LEGEND_POSITION = "Bottom";
    public static boolean MINI_HEATMAP = false;
    public static boolean MINI_HEATMAP_SINGLE = true;
    public static String SHOW_DATA_OPTION = "ALL";
    public static String RAW_PLOT_OPTION = "Plate View";

    //plate client
    public static boolean ADD_PLOT = false;
    public static boolean LEGO_VIEW = false;
    public static String PLATE_VIEW_OPTION = "single";
    public static String CONTROL_LABELS = "dbnisxyz";

    //correction
    public static String CORRECTION_ALGORITHM = "Default";
    public static String MULTI_VIEW_LAYOUT = "Continuous";
    public static boolean AUTO_CORRELATION_PLOT = false;
    public static int RANDOMIZATION_TIMES = 20;
    public static boolean SELECT_AND_MASK = true;
    public static int BLOCK_CORRECTION_OPTION = 0; //0=0 1=median 2=null

	//pivoting
	public static double ACTIVITY_DELTA = 0.0;
	public static double REPLICATE_PVALUE_CUTOFF = 0.05;
	public static boolean USE_CMPD_PLATE_ID2 = false;
	
	//auc
	public static double AUC_CONC_START = 1e-3; //um
	public static double AUC_CONC_END = 50; //um
	
	
	//query pre filter
	public static boolean QUERY_SELECT_CONVERT = true;
	public static boolean SAMPLE_ID_ROOT_MATCH = false;

}
