package com.openjfx;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.io.File;
import java.lang.Math;
import java.io.FileWriter;
import com.github.sh0nk.matplotlib4j.Plot;

public class Gascoyne_replication
{       
    static int DimensionalDepth = 10000;
    static Double Resolution = 0.002;
    static int NumberOfInputs = 10000;
    static Double SpatialUpperLimit = DimensionalDepth * Resolution;
    static Double InputTimeStep = 0.1;
    static Double TimeStep = 1.0;
    static int SampleRate = 1000;
    static int SampleResolution = 100;
    static Double FinalTime = InputTimeStep * 100000;
    static Double InputGaussianVariance = 35.0;
    static Double MemoryLossCoefficient = 0.0;
    public static void main(String[] args)
    {
        Double[] SpatialVectorField = VectorFieldInitialisation1Dimension(DimensionalDepth);
        Double[] PotentialVectorField = VectorFieldInitialisation1Dimension(DimensionalDepth);
        Function<Double,Double> EtaFunctionA = Eta -> 0.06 * Math.pow(Eta, 3) - (0.9 * Eta) - 0.9;
        Function<Double,Double> EtaFunctionC = Eta -> 0.08 * Math.pow(Eta, 3) +
                                                      40 * (Eta + 4)/Math.pow((Math.pow(Eta + 4,2) + 0.4),2) + 
                                                      50 * (Eta + 2)/Math.pow((Math.pow(Eta + 2,2) + 0.5),2) +
                                                      50 * (Eta - 0.5)/Math.pow((Math.pow(Eta - 0.5,2) + 1),2) +
                                                      300 * (Eta - 6)/Math.pow(Math.pow(Eta - 6,2) + 2, 2) +
                                                      300 * (Eta + 7)/Math.pow(Math.pow(Eta + 7,2) + 2, 2);
        Function<Double,Double> EtaFunctionD = Eta -> 0.2 * Eta + 6 * (Eta - 4) * Math.exp(-0.5 * Math.pow(Eta - 4,2)) +
                                                      3 * ((2 * Eta) - 10) * Math.exp(-Math.pow((Eta + 5),2)) +
                                                      8 * Eta * Math.exp(-Math.pow(Eta, 2)) +
                                                      5 * (8 * Eta - 56) * Math.exp(- 4 * Math.pow(Eta - 7, 2));
        Function<Double,Double> EtaFunctionE = Eta -> 0.2 * Eta + 5 * (0.4 * Eta + 0.4) * Math.exp(-0.2 * Math.pow(Eta + 1, 2)) +
                                                      7 * (6 * Eta + 42) * Math.exp(-3 * Math.pow(Eta + 7, 2)) +
                                                      6 * (Eta - 6) * Math.exp(-0.5 * Math.pow(Eta - 6, 2));

        Double[] DataSet = RevisedStimulusGeneration(NumberOfInputs,
                                             1.0,
                                             1.0,
                                             EtaFunctionA,
                                             InputTimeStep,
                                             FinalTime);
        File LocalOutputFile = new File("visualisation/src/main/resources/TotalSystemInput.csv");
        try
        {
            FileWriter writer = new FileWriter(LocalOutputFile);
            writer.append("V\n");
            for(int i = 0 ; i < DataSet.length; i++)
            {
                writer.append(DataSet[i] + "\n");
            }
            writer.close();
        }
        catch (Exception e)
        {   
            System.out.println("the following error has occured: " + e.toString());
        }
        Double minValue = Double.MAX_VALUE;
        for(Double element : DataSet)
        {
            if(element < minValue)
            
            {
                minValue = element;
            }
        }
        System.out.println(minValue);
        Double[] EvolvedSystem = PotentialEvolution(DataSet,
                                                    Resolution,
                                                    PotentialVectorField,
                                                    SpatialVectorField,
                                                    SampleRate,
                                                    SampleResolution);
        Double[] XData = Linspace(-DimensionalDepth/2*Resolution,DimensionalDepth/2*Resolution,DimensionalDepth);
        Double[] XDataDist = new Double[SpatialVectorField.length/SampleResolution];
        XDataDist = Linspace(minValue, -minValue, DimensionalDepth/SampleResolution);
        int[] BinnedDataset = RandomNumberGen.HistogramBinning(DataSet, SpatialVectorField.length/SampleResolution);
        Double[] DoubleBinnedData = new Double[SpatialVectorField.length/SampleResolution];
        for(int i = 0 ; i < XData.length; i++)
        {
 
        }
        File OutputFile = new File("visualisation/src/main/resources/SystemInput.csv");
        try
        {
            FileWriter writer = new FileWriter(OutputFile);
            writer.append("x, V\n");
            for(int i = 0 ; i < XData.length; i++)
            {
                if(i % SampleResolution == 0)
                {
                    DoubleBinnedData[i/SampleResolution] = -1.0*BinnedDataset[i/SampleResolution]/4000 ;
                    writer.append(XDataDist[i/SampleResolution] + "," + DoubleBinnedData[i/SampleResolution] + "\n");
                }  

            }
            writer.close();
        }
        catch (Exception e)
        {   
            System.out.println("the following error has occured: " + e.toString());
            for(int i = 0 ; i < XData.length; i++)
            {
                if(i % SampleResolution == 0)
                {
                    DoubleBinnedData[i/SampleResolution] = -1.0*BinnedDataset[i/SampleResolution];
                }                 
            }

        }
        LinePlot(XDataDist,
                 DoubleBinnedData,
                 "Evolved System",
                 "x",
                 "V",
                 "V(t,x)");                      
        LinePlot(XData,
                 EvolvedSystem,
                 "Evolved System",
                 "x",
                 "V",
                 "V(t,x)");
        LinePlot2Line(XData,
                      EvolvedSystem,
                      "Evolved system with input distribution",
                      "x",
                      "V(x,t)",
                      "V(x,t)",
                      XDataDist,
                      DoubleBinnedData,
                      "Distribution");
        System.out.println("MedianValue = " + EvolvedSystem[5000]);
        Double maxVal = 0.0;
        for(Double Data : EvolvedSystem)
        {
          if(Data < maxVal)
          {
            maxVal = Data;
          }
        }
        for(int i = 0; i < EvolvedSystem.length; i++)
        {
          if(EvolvedSystem[i] == maxVal)
          {
            System.out.println("central minima: " + i + "/" + EvolvedSystem.length);
          }  
        }
        System.out.println("ICs:");
        System.out.println("DimensionalDepth = " + DimensionalDepth + "\r\n" + //
                "Resolution = " + Resolution + "\r\n" + //
                "NumberOfInputs = " + NumberOfInputs+"\r\n" + //
                "SpatialUpperLimit = " + SpatialUpperLimit + "\r\n" + //
                "TimeStep = " + TimeStep);
    }
    public static Double[] GaussianFunction(Double Stimulus,
                                          Double Variance,
                                          Integer DataRange,
                                          int Dimensionality,
                                          double Increment)
    {
        Double[] DatapointList = new Double[DataRange];
        Double DataPoint = 1.0;
        Double SquareRootContents = (2* Math.PI * Math.pow(Variance, 2));
        Double GaussianPreFactor = 1/Math.pow(Math.sqrt(SquareRootContents),
                                              Dimensionality);
        for(int i = 0; i < DataRange; i++)
        {
            Double ModLimit = DataRange/2.0;
            DataPoint = GaussianPreFactor * Math.exp(-Math.pow(((i*Increment-(Increment*ModLimit)) - Stimulus), 2)
                                                        / (2*Math.pow(Variance, 2)));
            DatapointList[i] = -DataPoint;
        }
        return DatapointList;
    }
    public static Double[] gaussianStimulus(Double Stimulus,
                                        Double Variance,
                                        Integer Dimensionality,
                                        Double Increment)
    {
        File Datapoints = new File("Temp/gaussStimulus.txt");
        Double[] DatapointList = GaussianFunction(Stimulus, Variance, null, Dimensionality, Increment);
        try 
        {
            System.out.println(Datapoints.createNewFile());
            FileWriter DataFormatter = new FileWriter(Datapoints);
            DataFormatter.write("Dimensions :" + Dimensionality.toString() + "\n");
            
            for(int i = 0; i < DatapointList.length; i++)
            {
                for(Integer Dimension = 0; Dimension <= Dimensionality; Dimension ++)
                {
                    DataFormatter.write("Dimension: " + (Dimension + 1) +
                                        ", Location: " + i * Increment +
                                        ", Magnitude: " + DatapointList[i] +
                                        " ");
                }
                DataFormatter.write("\n");
            }
            DataFormatter.close();
        }
        catch (Exception e)
        {
            System.out.println("the following error has occured: " + e.toString());
        }
        return DatapointList;
    }

    public static Double[] PotentialVectorFieldChange(Integer DimensionalDepth,
                                                      Double Timestep,
                                                      Double[] StimulusVectorField,
                                                      Double MemoryLossCoefficient,
                                                      Double TotalTime,
                                                      Double[] PreviousPotentialVectorFieldIteration)
    {
        Double[] lPotentialVectorFieldChange = VectorFieldInitialisation1Dimension(DimensionalDepth);
        for(int i = 0; i < lPotentialVectorFieldChange.length; i++)
        {
            Double PreviousPotential = PreviousPotentialVectorFieldIteration[i];
            
            Double Datapoint = (1/TotalTime) * (PreviousPotential + StimulusVectorField[i])
                                - MemoryLossCoefficient * PreviousPotential;
            lPotentialVectorFieldChange[i] = Datapoint;
        }
        return lPotentialVectorFieldChange;
        
    }
    public static Double[] SpatialVectorFieldChange(Double RateOfLocalMinimaApproach,
                                                    Double[] PotentialField,
                                                    Double[] SystemNoise,
                                                    Double VectorFieldResolution)
    {
        Double[] SpatialPotentialDerivative1Dimension = VectorFieldInitialisation1Dimension(PotentialField.length);
        Double[] SpatialVectorFieldChange1Dimension = new Double[PotentialField.length];
        SpatialPotentialDerivative1Dimension[0] = (PotentialField[0]-PotentialField[1])/VectorFieldResolution;
        Double RHSDerivative = 0.0;
        Double LHSDerivative = 0.0;
        for(int i = 1; i < PotentialField.length-1; i++)
        {
            RHSDerivative = (PotentialField[i]-PotentialField[i+1])/VectorFieldResolution;
            LHSDerivative = (PotentialField[i-1]-PotentialField[i])/VectorFieldResolution;
            SpatialPotentialDerivative1Dimension[i] = (RHSDerivative + LHSDerivative)/2.0;
        }
        for(int i = 0; i < PotentialField.length; i++)
        {
            SpatialVectorFieldChange1Dimension[i] = -Math.sqrt(Math.pow(((RateOfLocalMinimaApproach * SpatialPotentialDerivative1Dimension[i]) + SystemNoise[i]),2));
        }
        return SpatialVectorFieldChange1Dimension;
    }
    public static Double[] VectorFieldInitialisation1Dimension(Integer DimensionalDepth)
    {
        Double[] VectorField1Dimension = new Double[DimensionalDepth];
        Arrays.fill(VectorField1Dimension, 0.0);
        return VectorField1Dimension;
    }
    public static Double[] PotentialEvolution(Double[] DataSet,
                                              Double Resolution,
                                              Double[] PotentialVectorField,
                                              Double[] SpatialVectorField,
                                              int SampleRate,
                                              int SampleResolution)
    {
        Double[] CurrentPotentialVectorField = VectorFieldInitialisation1Dimension(PotentialVectorField.length);
        Double[] CurrentSpatialVectorField = VectorFieldInitialisation1Dimension(SpatialVectorField.length);
        File OutputFile = new File("visualisation/src/main/resources/SystemLearning.csv");
        try
        {
            FileWriter writer = new FileWriter(OutputFile);
            writer.append("T, x, V(t-x)\n");
            for(int i = 0; i < DataSet.length; i++)
            {
                Double[] Potential = GaussianFunction(DataSet[i],
                                                      Resolution*InputGaussianVariance,
                                                      PotentialVectorField.length,
                                                      1,
                                                      Resolution);
                Double[] PotentialFieldChange = PotentialVectorFieldChange(PotentialVectorField.length, 
                                                                           1.0, 
                                                                           Potential, 
                                                                           MemoryLossCoefficient,
                                                                           1.0*DataSet.length, 
                                                                           CurrentPotentialVectorField);
                for(int j = 0; j < PotentialVectorField.length; j++)
                {
                    CurrentPotentialVectorField[j] += PotentialFieldChange[j];
                    if(i % SampleRate == 0)
                    {
                        if(j % SampleResolution == 0)
                        {
                        writer.append(i + "," + j*Resolution + "," + CurrentPotentialVectorField[j] + "\n");
                        }
                    }
                }
                Double[] SystemNoise = new Double[PotentialVectorField.length];
                Arrays.fill(SystemNoise, 0.0);
                Double[] SpatialFieldChange = SpatialVectorFieldChange(1.0,
                                                                       PotentialFieldChange,
                                                                       SystemNoise,
                                                                       Resolution);                                                      
                for(int k = 0; k < PotentialVectorField.length; k++)
                {
                    CurrentSpatialVectorField[k] += SpatialFieldChange[k];
                }
            }
            writer.close();
        }
        catch(Exception e)
        {
            System.out.println("The Following Error occurred: " + e.toString());
        }
        return CurrentPotentialVectorField;
    }
    public static void LinePlot(Double[] XDataPoints,
                                Double[] YDataPoints,
                                String Title,
                                String XAxisTitle,
                                String YAxisTitle,
                                String LineLabel)
    {
        List<Double> XDataList = Arrays.asList(XDataPoints);
        List<Double> YDataList = Arrays.asList(YDataPoints);
        Plot plt = Plot.create();
        plt.plot().add(XDataList, YDataList).label(LineLabel);
        plt.legend().loc("upper right");
        plt.title(Title);
        plt.xlabel(XAxisTitle);
        plt.ylabel(YAxisTitle);
        try
        {
            plt.show();
        }
        catch(Exception e)
        {

        }
    }
    public static void LinePlot2Line(Double[] XDataPoints,
                                     Double[] YDataPoints,
                                     String Title,
                                     String XAxisTitle,
                                     String YAxisTitle,
                                     String LineLabel,
                                     Double[] X2DataPoints,
                                     Double[] Y2DataPoints,
                                     String Line2Label)
    {
        List<Double> XDataList = Arrays.asList(XDataPoints);
        List<Double> YDataList = Arrays.asList(YDataPoints);
        List<Double> X2DataPointsList = Arrays.asList(X2DataPoints);
        List<Double> Y2DataPointsList = Arrays.asList(Y2DataPoints);
        Plot plt = Plot.create();
        plt.figure(Title);
        plt.plot().add(XDataList, YDataList).label(LineLabel);
        plt.plot().add(X2DataPointsList, Y2DataPointsList).label(Line2Label);
        plt.legend().loc("upper right");
        plt.title(Title);
        plt.xlabel(XAxisTitle);
        plt.ylabel(YAxisTitle);
        try
        {
            plt.show();
        }
        catch(Exception e)
        {

        }
    }
    public static Double[] Linspace(Double InitialValue,
                                    Double FinalValue, 
                                    int NumberOfDataPoints)
    {
        Double Interval = (FinalValue - InitialValue)/NumberOfDataPoints;
        Double[] LinspaceArray = new Double[NumberOfDataPoints];
        for(int i = 0; i < NumberOfDataPoints; i++)
        {
            LinspaceArray[i] = InitialValue + (i * Interval);
        }
        return LinspaceArray;
    }
    public static Double[] StimulusGeneration(int NumberOfGens,
                                              Double GaussianWhiteNoiseStrength,
                                              Double GaussianWhiteNoiseVariance,
                                              Function <Double, Double> StimulusFunction,
                                              Double TimeStep,
                                              Double FinalTime)
    {
        Random StimulusInputGenerator = new Random();
        Double[] StimulusList = new Double[NumberOfGens];
        Double[] GaussianWhiteNoise = RandomNumberGen.NumberGen(NumberOfGens, GaussianWhiteNoiseVariance);
        Double Stimulus = (StimulusInputGenerator.nextDouble() - 0.5) * SpatialUpperLimit * 1.2;
        Double Eta = Stimulus +  (-1.0 * StimulusFunction.apply(Stimulus) + GaussianWhiteNoiseStrength * GaussianWhiteNoise[0]) * TimeStep;
        StimulusList[0] = Eta;
        for(int i = 1; i < NumberOfGens; i++)
        {
            Stimulus = (StimulusInputGenerator.nextDouble() - 0.5)  * SpatialUpperLimit * 1.2;
            Eta = Stimulus - (StimulusFunction.apply(Stimulus) + GaussianWhiteNoiseStrength * GaussianWhiteNoise[i]) * TimeStep;
            for(int j = 0; j < FinalTime; j++)
            {
                Eta -= (StimulusFunction.apply(Eta) + GaussianWhiteNoiseStrength * GaussianWhiteNoise[i]) * TimeStep;
            }
            StimulusList[i] = Eta;
        }
        File Datapoints = new File("");
        File OutputFile = new File("visualisation/src/main/resources/SystemLearning.csv");
        return StimulusList;
    }
    public static Double[] RevisedStimulusGeneration(int NumberOfGens,
                                                     Double GaussianWhiteNoiseStrength,
                                                     Double GaussianWhiteNoiseVariance,
                                                     Function <Double, Double> StimulusFunction,
                                                     Double TimeStep,
                                                     Double Seed)
    {
        Random StimulusInputGenerator = new Random();
        Double[] StimulusList = new Double[NumberOfGens];
        Double[] GaussianWhiteNoise = RandomNumberGen.NumberGen(NumberOfGens, GaussianWhiteNoiseVariance);
        Double Stimulus = (StimulusInputGenerator.nextDouble() - 0.5) * SpatialUpperLimit * 1.2;
        Double Eta = Stimulus +  (-1.0 * StimulusFunction.apply(Stimulus) + GaussianWhiteNoiseStrength * GaussianWhiteNoise[0]) * TimeStep;
        StimulusList[0] = Eta;
        Double EtaDot = 0.0;
        for(int i = 1; i < NumberOfGens; i++)
        {
            EtaDot = -1.0 * StimulusFunction.apply(StimulusList[i-1]) + GaussianWhiteNoiseStrength * GaussianWhiteNoise[i] * TimeStep;
            StimulusList[i] = StimulusList[i-1] + EtaDot;
        }
        //File Datapoints = new File("");
        //File OutputFile = new File("visualisation/src/main/resources/SystemLearning.csv");
        return StimulusList;
    }
}