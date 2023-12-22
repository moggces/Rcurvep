package gov.nih.ncats.ifx.qhts.curvefitting;

import java.io.PrintStream;
import org.rosuda.JRI.RConsoleOutputStream;
import org.rosuda.JRI.Rengine;
import gov.nih.ncats.ifx.qhts.curvefitting.service.*;

public class CurveClassFit1 {

    private Rengine rEngine;

    public CurveClassFit1() {
        // Initialize Rengine and set System.out to RConsoleOutputStream
        rEngine = new Rengine();
        RConsoleOutputStream rConsoleOutputStream = new RConsoleOutputStream(rEngine, 0);
        PrintStream printStream = new PrintStream(rConsoleOutputStream);
        System.setOut(printStream);
    }

    public static void main(String[] args) {
        if (args.length != 4) {
            System.out.println("Usage: java CurveClassFit1 <conc> <resp> <classificationSD> <minYrange>");
            return;
        }

        CurveClassFit1 curveClassFit = new CurveClassFit1(); // Initialize CurveClassFit1

        // Parse command line arguments to double arrays and double values
        Double[] conc = parseDoubleArray(args[0]);
        Double[] resp = parseDoubleArray(args[1]);
        Double classificationSD = Double.parseDouble(args[2]);
        Double minYrange = Double.parseDouble(args[3]);

        // Perform curve fitting with provided arguments
        curveClassFit.performCurveFitting(conc, resp, classificationSD, minYrange);
    }

    // Method to perform curve fitting and output the result to redirected System.out
    private void performCurveFitting(Double[] conc, Double[] resp, Double classificationSD, Double minYrange) {
        CurveFittingService fittingService = new CurveFittingServiceImpl();
        FittingResult fittingResult = fittingService.doFit(conc, resp, classificationSD, minYrange);

        // Convert the fitting result to a String representation and print it
        String resultString = fittingResult.toString();
        System.out.println(resultString);
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

