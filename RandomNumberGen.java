package com.openjfx;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
public class RandomNumberGen 
{
    public static void main(String[] args)
    {
        
    }
    public static Double[] NumberGen(int NumberOfGens, Double Variance)
    {
        Double[] RandomNumberList = new Double[NumberOfGens];
        Double MinValue = Double.MAX_VALUE;
        Double MaxValue = Double.MIN_VALUE;
        for(int i = 0; i < NumberOfGens; i++)
        {
            RandomNumberList[i] = NormalisedNumberGen(Variance);
            if(RandomNumberList[i] < MinValue)
            {
                MinValue = RandomNumberList[i];
            }
            if(RandomNumberList[i] > MaxValue)
            {
                MaxValue = RandomNumberList[i];
            }
        }
        Double Midpoint = (MinValue-MaxValue)/2;
        for(int i = 0; i < NumberOfGens; i++)
        {
            //RandomNumberList[i] -= Midpoint;
        }
        return RandomNumberList;
    }
    public static Double NormalisedNumberGen(Double Variance)
    {
        Double nextNextGaussian = 0.0;
        Boolean haveNextNextGaussian = false;
        Random Generator  = new Random();
        // See Knuth, TAOCP, Vol. 2, 3rd edition, Section 3.4.1 Algorithm C.
        double v1, v2, s;
        do {
            v1 = 2 * Generator.nextDouble() - 1; // between -1 and 1
            v2 = 2 * Generator.nextDouble() - 1; // between -1 and 1
            s = v1 * v1 + v2 * v2;
        } while (s >= 1 || s == 0);
        double multiplier = Variance * StrictMath.sqrt(-2 * StrictMath.log(s)/s);
        nextNextGaussian = v2 * multiplier;
        haveNextNextGaussian = true;
        return v1 * multiplier;
    }
    public static int[] HistogramBinning(Double[] Dataset, int NumberOfBins)
    {
        Double MaxValue = Double.MIN_VALUE;
        Double MinValue = Double.MAX_VALUE;
        for(Double Data : Dataset)
        {
            if(Data > MaxValue)
            {
                MaxValue = Data;
            }
            if(Data < MinValue)
            {
                MinValue = Data;
            }
        }
        Double BinRange = (MaxValue - MinValue)/NumberOfBins;
        System.out.println("BinRange: " + BinRange + " maxval: " + MaxValue + " minval: " + MinValue);
        int[] BinnedData = new int[NumberOfBins+1];
        int count = 0;
        for(double i = MinValue; i <= MaxValue ; i += BinRange)
        {
            for(int j = 0; j < Dataset.length; j++)
            {
                if(Dataset[j] < (i + BinRange) &&
                   Dataset[j] > i)
                {
                    BinnedData[count] += 1;
                }
            }
            count += 1;
        }
        return BinnedData;
    }
    public static Double[]OrsteinUhlenbeckPathGeneration(int NumberOfGens,
                                                         Double Theta,
                                                         Double Mu,
                                                         Double Sigma,
                                                         Double InitialValue,
                                                         Double InitialTime,
                                                         Double FinalTime)
    {
        Double[] OutputNumbers = new Double[NumberOfGens];
        Double DeltaT = (FinalTime - InitialTime)/NumberOfGens;
        Double RootDeltaT = Math.sqrt(DeltaT);
        Double DeltaW = LocalnextGaussian(RootDeltaT);
        Double FirstRandom = InitialValue  + (Theta * (Mu - InitialValue) * DeltaT) + (Sigma  * DeltaW);
        OutputNumbers[0] = FirstRandom;
        for(int i = 1; i < NumberOfGens; i++)
        {
            DeltaW = LocalnextGaussian(RootDeltaT);
            OutputNumbers[i] = OutputNumbers[i-1] +  (Theta * (Mu - OutputNumbers[i-1]) * DeltaT) + (Sigma * DeltaW);;
        }
        return OutputNumbers;
    }
    public static Double[]NoiseDistribution(int NumberOfDatapoints,
                                            Double Theta,
                                            Double Mu,
                                            Double Sigma,
                                            Double InitialValue,
                                            Double InitialTime,
                                            Double FinalTime,
                                            Double MinValue,
                                            Double MaxValue)
    {
        double DeltaT = FinalTime - InitialTime;
        if(DeltaT % 2 < 0.5)
        {
            DeltaT += 0.5;
        }
        int DelatTint = (int)DeltaT;
        Double Increment = (MaxValue - MinValue)/NumberOfDatapoints;
        Double[] ProbabilityDistribution = new Double[NumberOfDatapoints];
        Double D = Math.pow(Sigma, 2)/2;
        Double RootDenominator = 2 * Math.PI * D * (1 - Math.exp(-2 * Theta * (DeltaT)));
        Double DistributionPreFactor = Math.sqrt(Theta/RootDenominator);
        Double ExponentialDenominatorConstant = - (Theta / 2 * D) * (1 /(1 - (Math.exp(-2 * Theta * DeltaT))));
        Double ExponentialNumeratorConstant = Math.exp(- Theta * DelatTint);
        for(int i = 0; i < NumberOfDatapoints; i++)
        {
            ProbabilityDistribution[i] = DistributionPreFactor * Math.exp(ExponentialDenominatorConstant * Math.pow(((i*Increment) - (InitialValue * ExponentialNumeratorConstant)),2));
        }
        return ProbabilityDistribution;
    }
    public static Double[]NoiseWiener(int NumberOfGens,
                                      Double Theta,
                                      Double Sigma,
                                      Double InitialValue,
                                      Double InitialTime,
                                      Double FinalTime)
    {
        Double[] WienerNoise = new Double[NumberOfGens];
        Double CurrentX = 0.0;
        Double WienerS  = 0.0;
        Double Timestep = (FinalTime - InitialTime)/NumberOfGens;
        Double CurrentTime = 0.0;
        for(int i = 0; i < NumberOfGens; i++)
        {
            CurrentTime = i * Timestep;
            WienerS = WienerProcess(Theta, Sigma, CurrentTime, InitialValue);
            CurrentX = (Sigma/Math.sqrt(2 * Theta)) * Math.exp(-Theta * CurrentTime) * WienerS;
            WienerNoise[i] = CurrentX;
        }
        return WienerNoise;

    }
    public static Double WienerProcess(Double Theta,Double Sigma, Double Time, Double InitialValue)
    {
        Double SigmaTheta = Math.sqrt(2 * Theta)/ Sigma;
        Double ThetaExponentTerm = Math.exp(2* Theta * Time);
        Double WienerS = SigmaTheta * Math.sqrt(ThetaExponentTerm) * InitialValue * Math.log(ThetaExponentTerm)
                        /2 * Theta;
        return WienerS;  
    }
    private double nextNextGaussian;
    private boolean haveNextNextGaussian = false;
    public static double LocalnextGaussian(Double Variance)
    {
        double nextNextGaussian = 0.0;
        boolean haveNextNextGaussian = false;
        Random Generator  = new Random();
        // See Knuth, TAOCP, Vol. 2, 3rd edition, Section 3.4.1 Algorithm C.
        if (haveNextNextGaussian) {
            haveNextNextGaussian = false;
            return nextNextGaussian;
        } else {
            double v1, v2, s;
            do {
                v1 = 2 * Generator.nextDouble() - 1; // between -1 and 1
                v2 = 2 * Generator.nextDouble() - 1; // between -1 and 1
                s = v1 * v1 + v2 * v2;
            } while (s >= 1 || s == 0);
            double multiplier = Variance * StrictMath.sqrt(-2 * StrictMath.log(s)/s);
            nextNextGaussian = v2 * multiplier;
            haveNextNextGaussian = true;
            return v1 * multiplier;
        }
    }
}
