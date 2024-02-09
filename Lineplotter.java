package com.openjfx;

import java.io.BufferedReader;
import java.io.FileReader;
import javafx.application.Application;
import javafx.scene.Group;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.stage.Stage;

public class Lineplotter extends Application{
   @Override     
   public void start(Stage primaryStage) throws Exception
   {
      //Defining X axis  

      XYChart.Series series = new XYChart.Series(); 
      series.setName("RandomNumbers");
      String FilePath ="visualisation/src/main/resources/Histogram.txt";
      BufferedReader reader = new BufferedReader(new FileReader(FilePath));
      String currentLine = reader.readLine();
      String[] AxisTitles = currentLine.split(",");
      Double Xlim = Double.MIN_VALUE;
      Double Ylim = Double.MIN_VALUE;
      Double XlimMin = Double.MAX_VALUE;
      Double YlimMin = Double.MAX_VALUE;
      currentLine = reader.readLine();
      while(currentLine != null)
      {
         String[] StringData = currentLine.split(",");
         Double[] Data = {Double.parseDouble(StringData[0]), Double.parseDouble(StringData[1])};
         series.getData().add(new XYChart.Data(Data[1], Data[0]));
         currentLine = reader.readLine();
         if(Data[1] > Xlim)
         {
            Xlim = Data[1];
         }
         if(Data[0] > Ylim)
         {
            Ylim = Data[0];
         }
         if(Data[1] < XlimMin)
         {
            XlimMin = Data[1];
         }
         if(Data[0] < YlimMin)
         {
            YlimMin = Data[0];
         }
      }

      reader.close();      
      NumberAxis xAxis = new NumberAxis(XlimMin, Xlim, Xlim/10); 
      xAxis.setLabel(AxisTitles[1]);
      //Defining y axis 
      NumberAxis yAxis = new NumberAxis(0,1.5*Ylim , Ylim/10); 
      yAxis.setLabel(AxisTitles[0]);
      LineChart linechart = new LineChart(xAxis, yAxis);
      //Setting the data to Line chart    
      linechart.getData().add(series);
      Group root = new Group(linechart);
      Scene scene = new Scene(root ,1200, 600);
      primaryStage.setTitle("Random NumberFrequency");
      primaryStage.setScene(scene);
      primaryStage.show();
   }    
   public static void main(String args[]){   
    launch(args);   
   }
   public void LinePlot(Stage primaryStage,
                        String PlotTitle,
                        String XAxisLabel,
                        String YAxisLabel,
                        Double[] XDataPoints,
                        Double[] YDataPoints) throws Exception
   {
      //Defining X axis 
      XYChart.Series series = new XYChart.Series(); 
      series.setName("RandomNumbers");
      Double Xlim = Double.MIN_VALUE;
      Double Ylim = Double.MIN_VALUE;
      Double XlimMin = Double.MAX_VALUE;
      Double YlimMin = Double.MAX_VALUE;
      for(int i =0; i< XDataPoints.length; i++)
      {
         series.getData().add(new XYChart.Data(XDataPoints[i], YDataPoints[0]));
         if(XDataPoints[i] > Xlim)
         {
            Xlim = XDataPoints[i];
         }
         if(YDataPoints[i] > Ylim)
         {
            Ylim = YDataPoints[i];
         }
         if(XDataPoints[i] < XlimMin)
         {
            XlimMin = XDataPoints[i];
         }
         if(YDataPoints[i] < YlimMin)
         {
            YlimMin = YDataPoints[i];
         }
      }   
      NumberAxis xAxis = new NumberAxis(XlimMin, Xlim, Xlim/10); 
      xAxis.setLabel(XAxisLabel);
      //Defining y axis 
      NumberAxis yAxis = new NumberAxis(0,1.5*Ylim , Ylim/10); 
      yAxis.setLabel(YAxisLabel);
      LineChart linechart = new LineChart(xAxis, yAxis);
      //Setting the data to Line chart    
      linechart.getData().add(series);
      Group root = new Group(linechart);
      Scene scene = new Scene(root ,1200, 600);
      primaryStage.setTitle(PlotTitle);
      primaryStage.setScene(scene);
      primaryStage.show();
   } 
}
