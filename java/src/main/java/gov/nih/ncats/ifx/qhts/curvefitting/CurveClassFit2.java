package gov.nih.ncats.ifx.qhts.curvefitting;

import gov.nih.ncats.ifx.qhts.curvefitting.service.*;

public class CurveClassFit2 {

    public static void main(String[] args) {
        if (args.length != 4) {
            System.out.println("Usage: java CurveClassFit2 <conc> <resp> <classificationSD> <minYrange>");
            return;
        }

        // Parse command line arguments to double arrays and double values
        Double[] conc = parseDoubleArray(args[0]);
        Double[] resp = parseDoubleArray(args[1]);
        Double classificationSD = Double.parseDouble(args[2]);
        Double minYrange = Double.parseDouble(args[3]);

        // Perform curve fitting with provided arguments and print the result
        String result = performCurveFitting(conc, resp, classificationSD, minYrange);
        System.out.println(result);
    }

    // Method to perform curve fitting and return the result as a String
    public static String performCurveFitting(Double[] conc, Double[] resp, Double classificationSD, Double minYrange) {
        CurveFittingService fittingService = new CurveFittingServiceImpl();
        FittingResult fittingResult = fittingService.doFit(conc, resp, classificationSD, minYrange);
        
        // Convert the fitting result to a String representation
        return fittingResult.toString(); // Assuming FittingResult has overridden toString() method
    }

    // Method to parse space-separated string of doubles into a Double array
    private static Double[] parseDoubleArray(String arg) {
        String[] values = arg.split("\\s+");
        Double[] result = new Double[values.length];
        for (int i = 0; i < values.length; i++) {
            result[i] = Double.parseDouble(values[i]);
        }
        return result;
    }
}
