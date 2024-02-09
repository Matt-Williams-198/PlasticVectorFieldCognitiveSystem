package com.openjfx;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.lang.Math;
import java.io.FileWriter;
import java.io.OutputStream;

import com.github.sh0nk.matplotlib4j.Plot;
import com.openjfx.Gascoyne_replication;    
import java.lang.reflect.*;
import java.net.URI;
public class MultiDimensionalSystem
{
    public static void main(String[] args)
    {
        Object abc = nArray(3,10);
        System.out.println(abc);
    }
    public static List<Double[]> MultiDimensionalVectorFieldInitialisation(int NumberOfDims, int DimensionalDepth)
    {
        Double[] ContentArray = Gascoyne_replication.VectorFieldInitialisation1Dimension(DimensionalDepth);
        //List<Double[]> MultiDimensionalVectorField = new ArrayList<Double[]>(DimensionalDepth);
        ArrayList<ArrayList<Double>> MultiDimensionalVectorField = new ArrayList<>();
        //Array.set
        return new ArrayList();
    }

    public static Object nArray(int n, int DimensionalDepth)
    {
        int[]dim = new int [DimensionalDepth];
        Arrays.fill(dim, DimensionalDepth);
        return Array.newInstance(int.class,dim);
    }
}
