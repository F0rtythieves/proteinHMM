import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.Map;

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

    public static int[] add(int[] first, int[] second) {
        int length = first.length < second.length ? first.length : second.length;
        int[] result = new int[length];
        for (int i = 0; i < length; i++) {
            result[i] = first[i] + second[i]; 
        } 
        return result; 
    }

    public static int[] countFragments(String stateSeq){
        int phobic = stateSeq.split("O", -1).length-1;
        int philic = stateSeq.split("I", -1).length-1;
        int mixed = stateSeq.split("M", -1).length-1;
        return new int[] {phobic, philic, mixed};
    }

    public static Map countLengths(ArrayList<String> annotations, char region){
        Map<Integer, Integer> dist = new HashMap<Integer, Integer>();
        //[hydrophobic, hydrophilic, mixed]

        for (int i = 0; i < annotations.size(); i++){
            String s = annotations.get(i);
            int length = 0;
            for (int c = 0; c < s.length(); c++){
                if (s.charAt(c) == region){
                    length++;
                }
                else{
                    if (length != 0){
                        if (dist.containsKey(length)){
                            dist.put(length, dist.get(length) + 1);
                        }
                        else{
                            dist.put(length, 1);
                        }
                        length = 0;
                    }
                }
            }
        }
        return dist;
    }

    public static Map countFrequency(String sequence){
        Map<Character, Integer> frequencies = new HashMap<Character, Integer>();

        for (int i = 0; i < sequence.length(); i++){
            char c = sequence.charAt(i);
            if (frequencies.containsKey(c)){
                frequencies.put(c, frequencies.get(c) + 1);
            }
            else{
                frequencies.put(c, 1);
            }
        }

        return frequencies;
    }

    public static void main(String args[]) throws IOException {
        System.out.println("O - Hydrophobic region");
        System.out.println("I - Hydrophilic region");
        System.out.println("M - Mixed region");

        BufferedReader in = new BufferedReader(new FileReader("hw3_proteins.fa"));
        ArrayList<String> stateSequences = new ArrayList<String>();
        ArrayList<String> proteinSequences = new ArrayList<String>();

        String line;
        StringBuilder proteinSeq = new StringBuilder();
        Boolean newProtein = false;
        while ((line=in.readLine()) != null){
            if (line.startsWith(">")){
                // System.out.println(line);

                if (newProtein){
                    String stateSeq = getStateSeq(proteinSeq.toString());
                    // System.out.println(stateSeq);
                    stateSequences.add(stateSeq);
                    proteinSequences.add(proteinSeq.toString());
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
        stateSequences.add(stateSeq);
        proteinSequences.add(proteinSeq.toString());


        // Map oDist = countLengths(stateSequences, 'O');
        // Map iDist = countLengths(stateSequences, 'I');
        // Map mDist = countLengths(stateSequences, 'M');

        StringBuilder oAminoAcids = new StringBuilder();
        StringBuilder iAminoAcids = new StringBuilder();
        StringBuilder mAminoAcids = new StringBuilder();

        for (int i = 0; i < proteinSequences.size(); i++){
            String sSeq = stateSequences.get(i);
            String pSeq = proteinSequences.get(i);
            for (int j = 0; j < sSeq.length(); j++){
                if (sSeq.charAt(j) == 'O'){
                    oAminoAcids.append(pSeq.charAt(j));
                }
                else if (sSeq.charAt(j) == 'I'){
                    iAminoAcids.append(pSeq.charAt(j));
                }
                else if (sSeq.charAt(j) == 'M'){
                    mAminoAcids.append(pSeq.charAt(j));
                }
            }

        }

        Map oFrequencies = countFrequency(oAminoAcids.toString());
        Map iFrequencies = countFrequency(iAminoAcids.toString());
        Map mFrequencies = countFrequency(mAminoAcids.toString());

        System.out.println(oAminoAcids.length());
        System.out.println(oFrequencies);

        System.out.println(iAminoAcids.length());
        System.out.println(iFrequencies);

        System.out.println(mAminoAcids.length());
        System.out.println(mFrequencies);

    }
}