import java.io.*;
import java.util.*;
import java.util.HashMap;
import java.util.Map;

public class MotifDiscovery {

    public static void printMatrix(int[][] m){
        for (int i =0; i < m.length; i++){
            for (int j = 0; j < m[i].length; j++){
                System.out.print(m[i][j] + "|");
            }
            System.out.println("\n_");
        }
        System.out.println("************************");
    }

    public static ArrayList<String> loadFa(String filePath) throws IOException{
        BufferedReader gata = new BufferedReader(new FileReader(filePath));
        ArrayList<String> gataSequences = new ArrayList<String>();

        String line;
        StringBuilder gataSeq = new StringBuilder();
        Boolean newProtein = false;
        while ((line=gata.readLine()) != null){
            if (line.startsWith(">")){
                if (newProtein){
                    gataSequences.add(gataSeq.toString());
                    newProtein = false;
                    gataSeq = new StringBuilder();
                }
                newProtein = true;
            }
            else{
                gataSeq.append(line);
            }   
        }
        gataSequences.add(gataSeq.toString());

        return gataSequences;
    }

    public static int countMatches(String consensus, ArrayList<String> sequences){
        int matches = 0;
        for (int i = 0; i < sequences.size(); i++){
            // System.out.println("Processing sequence: " + (i+1));
            String seq = sequences.get(i);
            for (int j = 0; j < seq.length() - (consensus.length() - 1); j++){
                int tempMatches = 0;
                for (int k = 0; k < consensus.length(); k++){
                    int seqPos = k+j;
                    if (consensus.charAt(k) == seq.charAt(seqPos)){
                        tempMatches++;
                    }
                    else if (consensus.charAt(k) == 'u'){
                        if (seq.charAt(seqPos) == 'A' || seq.charAt(seqPos) == 'G'){
                            tempMatches++;
                        }
                    }
                    else if (consensus.charAt(k) == 'y'){
                        if (seq.charAt(seqPos) == 'C' || seq.charAt(seqPos) == 'T'){
                            tempMatches++;
                        }
                    }
                    else if (consensus.charAt(k) == 'n'){
                        tempMatches++;
                    }
                    else{
                        break;
                    }
                }
                if (tempMatches == consensus.length()){
                    matches++;
                }
            }
        }
        return matches;
    }

    public static void main(String args[]) throws IOException{
        ArrayList<String> gataSequences = loadFa("GATA2_chr1.fa");
        ArrayList<String> notGataSequences = loadFa("not_GATA2_chr1.fa");

        // u - A|G
        // y - C|T
        // n - A|C|G|T
        char[] alphabet = new char[]{'A', 'C', 'G', 'T', 'u', 'y', 'n'};

        int consensusLength = 6;

        StringBuilder cSequence = new StringBuilder("      ");
        ArrayList<String> consensusSequences = new ArrayList<String>();

        int[] pos = new int[consensusLength];
        int totalSequences = (int) Math.pow(consensusLength, alphabet.length);

        for (int i =0; i < totalSequences; i++){
            for (int j = 0; j < consensusLength; j++){
                if (pos[j] == alphabet.length){
                    pos[j] = 0;
                    if (j + 1 < consensusLength){
                        pos[j+1]++;
                    }
                }
                cSequence.setCharAt(j, alphabet[pos[j]]);
            }
            pos[0]++;
            consensusSequences.add(cSequence.toString());
        }

        double bestZscore = -1000000000.0;
        String bestConsensus = "";
        for (int i = 0; i < consensusSequences.size(); i++){
            long startTime = System.currentTimeMillis();

            String cSeq = consensusSequences.get(i);
            int consensusMatches = countMatches(cSeq, gataSequences);
            int expectedMatches = countMatches(cSeq, notGataSequences);

            double zscore = ((double)(consensusMatches - expectedMatches))/Math.sqrt(expectedMatches);
            if (zscore > bestZscore){
                bestZscore = zscore;
                bestConsensus = cSeq;
            }
            if ((i+1) % 100 == 0){
                System.out.print(i + "/" + (consensusSequences.size()-1) + ": ");
                System.out.println(bestZscore);
                System.out.println(bestConsensus);
                double elapsedTime = (System.currentTimeMillis() - startTime)/1000.0;
                double timeLeft = elapsedTime * (consensusSequences.size() - i -1);
                System.out.println(timeLeft/3600 + " hours left.");
            }
            
        }


        // System.out.println(consensusSequences.size());

        // System.out.println(gataSequences.size());
        // System.out.println(notGataSequences.size());

    }
}