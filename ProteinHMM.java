import java.io.*;
import java.util.*;

public class ProteinHMM {
    public static double[] initialProb = new double[]{
        1/3.0,
        1/3.0,
        1/3.0};
    public static double[][] emissionProb = new double[][]{
        {0.6/8, 0.4/12},
        {0.2/8, 0.8/12},
        {1.0/20, 1.0/20}};
    public static double[][] transitionProb = new double[][]{
        {4.0/5, 0.2/5, 0.8/5},
        {0.3/8, 7.0/8, 0.7/8},
        {0.5/7, 0.5/7, 6.0/7}};

    public static void printMatrix(int[][] m){
        for (int i =0; i < m.length; i++){
            for (int j = 0; j < m[i].length; j++){
                System.out.print(m[i][j] + "|");
            }
            System.out.println("\n_");
        }
        System.out.println("************************");
    }

    public static int getIndex(char a){
        if ("AVILMFYW".indexOf(a) != -1){
            return 0;
        }
        else if ("RHKDESTNQCGP".indexOf(a) != -1){
            return 1;
        }
        return -1;
    }

    public static String getStateSeq(String protein){
        double[][] v = new double[protein.length()][3];
        int[][] vMaxState = new int[protein.length()][3];

        v[0][0] = initialProb[0] * emissionProb[0][getIndex(protein.charAt(0))];
        v[0][1] = initialProb[1] * emissionProb[1][getIndex(protein.charAt(0))];
        v[0][2] = initialProb[2] * emissionProb[2][getIndex(protein.charAt(0))];

        vMaxState[0][0] = 0;
        vMaxState[0][1] = 1;
        vMaxState[0][2] = 2;

        for (int i = 1; i < protein.length(); i++){
            char aa = protein.charAt(i);

            for (int s = 0; s < 3; s++){
                int maxState = 0;
                double prevStatProb = 0.0;
                double ho = v[i-1][0] * transitionProb[0][s]; // hydrophobic
                // System.out.println(ho);
                double hi = v[i-1][1] * transitionProb[1][s]; // hydrophilic
                // System.out.println(hi);
                double m = v[i-1][2] * transitionProb[2][s];  // mixed
                // System.out.println(m);

                if ((ho > hi) && (ho > m)){
                    maxState = 0;
                    prevStatProb = ho;
                }
                else if ((hi > ho) && (hi > m)){
                    maxState = 1;
                    prevStatProb = hi;
                }
                else{
                    maxState = 2;
                    prevStatProb = m;
                }

                // System.out.println(maxState);

                vMaxState[i][s] = maxState;

                v[i][s] = prevStatProb * emissionProb[s][getIndex(aa)];
            }            
        }
        // printMatrix(v);
        // printMatrix(vMaxState);

        int maxIndex = 0;
        double temp = v[protein.length() - 1][0];
        for (int i = 1; i < 3; i++){
            if (v[protein.length() - 1][i] > temp){
                temp = v[protein.length() - 1][i];
                maxIndex = i;
            }
        }

        String states = "OIM";
        StringBuilder stateSequence = new StringBuilder();
        // System.out.println(maxIndex);
        int tempIndex = maxIndex;
        for (int i = protein.length()-1; i >=0; i--){
            stateSequence.append(states.charAt(vMaxState[i][tempIndex]));
            tempIndex = vMaxState[i][tempIndex];
        }
        return stateSequence.reverse().toString();
    }

    public static void main(String args[]) throws IOException {
        System.out.println("O - Hydrophobic region");
        System.out.println("I - Hydrophilic region");
        System.out.println("M - Mixed region");

        BufferedReader in = new BufferedReader(new FileReader("hw3_proteins.fa"));

        String line;
        StringBuilder proteinSeq = new StringBuilder();
        int i = 0; 
        Boolean newProtein = false;
        while ((line=in.readLine()) != null){
            if (line.startsWith(">")){
                System.out.println(line);

                if (newProtein){
                    String stateSeq = getStateSeq(proteinSeq.toString());
                    System.out.println(stateSeq);
                    // System.out.println(proteinSeq);
                    newProtein = false;
                    proteinSeq = new StringBuilder();
                }
                newProtein = true;
            }
            else{
                proteinSeq.append(line);
            }   
        }
        // To account for the last one in file.
        String stateSeq = getStateSeq(proteinSeq.toString());
        System.out.println(stateSeq);

    }
}