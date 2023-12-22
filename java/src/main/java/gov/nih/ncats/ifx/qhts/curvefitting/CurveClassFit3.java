package gov.nih.ncats.ifx.qhts.curvefitting;

import gov.nih.ncats.ifx.qhts.curvefitting.service.*;

public class CurveClassFit3 {

    public static void main(String[] args) {
        if (args.length != 4) {
            System.out.println("Usage: java CurveClassFit3 <conc> <resp> <classificationSD> <minYrange>");
            return;
        }

        // Parse command line arguments to double arrays and double values
        double[] conc = parseDoubleArray(args[0]);
        double[] resp = parseDoubleArray(args[1]);
        double classificationSD = Double.parseDouble(args[2]);
        double minYrange = Double.parseDouble(args[3]);

        // Perform curve fitting with provided arguments and print the result
        String result = performCurveFitting(conc, resp, classificationSD, minYrange);
        System.out.println(result);
    }

    // Method to perform curve fitting and return the result as a String
    public static String performCurveFitting(double[] conc, double[] resp, double classificationSD, double minYrange) {
        CurveFittingService fittingService = new CurveFittingServiceImpl();
        FittingResult fittingResult = fittingService.doFit(toDoubleArray(conc), toDoubleArray(resp), classificationSD, minYrange);
        
        // Convert the fitting result to a String representation
        return fittingResult.toString(); // Assuming FittingResult has overridden toString() method
    }

    // Method to parse space-separated string of doubles into a double array
    private static double[] parseDoubleArray(String arg) {
        String[] values = arg.split("\\s+");
        double[] result = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = Double.parseDouble(values[i]);
        }
        return result;
    }

    // Convert double[] to Double[]
	private static Double[] toDoubleArray(double[] arr) {
		Double[] result = new Double[arr.length];
		for (int i = 0; i < arr.length; i++) {
			result[i] = arr[i]; // Autoboxing: converting double to Double
		}
		return result;
	}
}
