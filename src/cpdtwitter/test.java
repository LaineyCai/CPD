package cpdtwitter;

import java.io.File;
import java.io.IOException;

import de.bwaldvogel.liblinear.Feature;
import de.bwaldvogel.liblinear.FeatureNode;
import de.bwaldvogel.liblinear.Linear;
import de.bwaldvogel.liblinear.Parameter;
import de.bwaldvogel.liblinear.Problem;
import de.bwaldvogel.liblinear.SolverType;

public class test {
	public static void main(String args[]) throws IOException
	{
	    //loading train data
        Feature[][] featureMatrix = new Feature[5][];
        double[][] dv = new double[][]{{0,0.1,0.2,0,0},{0,0.1,0.3,-1.2,0},{0.4,0,0,0,0},{0,0.1,0,1.4,0.5},{-0.1,-0.2,0.1,-1.1,0.1}};
//        Feature[] featureMatrix1 = { new FeatureNode(2, 0.1), new FeatureNode(3, 0.2) };
//        Feature[] featureMatrix2 = { new FeatureNode(2, 0.1), new FeatureNode(3, 0.3), new FeatureNode(4, -1.2)};
//        Feature[] featureMatrix3 = { new FeatureNode(1, 0.4) };
//        Feature[] featureMatrix4 = { new FeatureNode(2, 0.1), new FeatureNode(4, 1.4), new FeatureNode(5, 0.5) };
//        Feature[] featureMatrix5 = { new FeatureNode(1, -0.1), new FeatureNode(2, -0.2), new FeatureNode(3, 0.1), new FeatureNode(4, -1.1), new FeatureNode(5, 0.1) };
//        featureMatrix[0] = featureMatrix1;
//        featureMatrix[1] = featureMatrix2;
//        featureMatrix[2] = featureMatrix3;
//        featureMatrix[3] = featureMatrix4;
//        featureMatrix[4] = featureMatrix5;
        int M = dv.length, N = dv[0].length;
        for(int i=0; i<M; i++)
        {
        	featureMatrix[i] = new Feature[N];
        	for(int j=0; j<N; j++)
        	{
        		featureMatrix[i][j]=new FeatureNode(j+1,dv[i][j]);
        	}
        }
        
        //loading target value
        double[] targetValue = {1,-1,1,-1,1};
        
        Problem problem = new Problem();
        problem.l = 5; // number of training examples：训练样本数
        problem.n = 5; // number of features：特征维数
        problem.x = featureMatrix; // feature nodes：特征数据
        problem.y = targetValue; // target values：类别

        SolverType solver = SolverType.L2R_LR; // -s 0
        double C = 1.0;    // cost of constraints violation
        double eps = 0.01; // stopping criteria
            
        Parameter parameter = new Parameter(solver, C, eps);
        de.bwaldvogel.liblinear.Model model = Linear.train(problem, parameter);
        File modelFile = new File("model");
        model.save(modelFile);
        // load model or use it directly
        model = de.bwaldvogel.liblinear.Model.load(modelFile);
        
        double[] w = model.getFeatureWeights();
        for(double dw:w)
        {
        	System.out.println(dw);
        }

        Feature[] testNode = { new FeatureNode(1, 0.4), new FeatureNode(3, 0.3) };//test node
        double prediction = Linear.predict(model, testNode);
        System.out.print("classification result: "+prediction);
	}

}
