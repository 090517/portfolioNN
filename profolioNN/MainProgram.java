package profolioNN;

import java.io.File;
import java.io.IOException;

import activationFunctions.ActivationFunction;
import activationFunctions.LeakyRelu;
import activationFunctions.Relu;
import activationFunctions.Sigmoid;
import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import neuralNet.MLNN;

public class MainProgram {
	public static void main(String[] arg) {
		long startTime = System.currentTimeMillis();
		System.out.println("start");

		double[][][] inputData = new double[250][60][14];
		double[][][] charData = new double[250][11][10];
		double[][][] idealData = new double[250][12][12];
		double[][] idealWeights = new double[250][10];
		double[] riskFree = new double[250];
		
		// read in ALL
		try {
			idealWeights = read0("idealWeights.xls");
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		String inputFile;
		for (int dataFile = 1; dataFile <= 250; dataFile++) {
			inputFile = new String(dataFile + "_data_practice.xls");
			try {
				inputData[dataFile - 1] = read1(inputFile);
				charData[dataFile - 1] = read2(inputFile);
				idealData[dataFile - 1] = read3(inputFile);
				riskFree[dataFile-1]= read4(inputFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		//double[] equalWeights= {.1,.1,.1,.1,.1,.1,.1,.1,.1,.1};
		//System.out.println(sharpCalc(equalWeights, idealData[0], 0.0071));
		//double[] unequalWeights= {-.1,.2,-.1,.2,-.1,.2,-.1,.2,-.1,.2};
		//System.out.println(sharpCalc(unequalWeights, idealData[0], 0.0071));
		
		double averageSharp=0;
		for (int i=0; i<150; i++) {
			averageSharp+=sharpCalc(idealWeights[i], idealData[i], riskFree[i]);
		}
		System.out.println("Average training data max sharp is:" + averageSharp/150);
		
		averageSharp=0.0;
		
		for (int i=150; i<250; i++) {
			averageSharp+=sharpCalc(idealWeights[i], idealData[i], riskFree[i]);
		}
		System.out.println("Average validation data max sharp is:" + averageSharp/100);
		
		double[][] inputArray = new double[250][950];
		for (int i = 0; i < 250; i++) {
			for (int j = 0; j < 60; j++) {
				for (int k = 0; k < 14; k++) {
					inputArray[i][j * 14 + k] = inputData[i][j][k];
					//System.out.print(inputData[i][j][k] + "\t");
				}
				 //System.out.println();
			}
			for (int j = 0; j < 11; j++) {
				for (int k = 0; k < 10; k++) {
					inputArray[i][840 + j * 10 + k] = charData[i][j][k];
					//System.out.print(charData[i][j][k] + "\t");
				}
				//System.out.println();
			}
		}
		for (int i=0; i<950; i++) {
			//System.out.print(inputArray[0][i] + "\t");
			if((i+1)%14==0&&i<840) {
				//System.out.println();
			}
			if((i+1)%10==0&&i>=840) {
				//System.out.println();
			}
		}
		
		
		//print3d(inputData);
		//print3d(charData);
		//print3d(idealData);

		// create neural net	
		
		int hiddenlayers=3;
		int[] hiddenArray = new int[hiddenlayers];
		hiddenArray[0] = 950;
		hiddenArray[1] = 950;
		hiddenArray[2] = 950;
		
		double leakyReluRate=.001;
		ActivationFunction[] AF = new ActivationFunction[hiddenlayers+1];
		for (int i = 0; i < AF.length; i++) {
			
			AF[i] = new Sigmoid();
			//AF[i] = new LeakyReLU(leakyReluRate);
		}
		AF[hiddenlayers] = new LeakyRelu(leakyReluRate);

		boolean[] bias = new boolean[hiddenlayers+1];
		for (int i = 0; i < AF.length; i++) {
			bias[i] = true;
		}
		
		double learningRate = .7;
		int trainingRounds=100;
		
		MLNN portfolioNN = new MLNN(hiddenArray, 950, 10, AF, bias, learningRate);
		// portfolioNN=MLNN.readNN("portfolioNNwithValdiation-A" + portfolioNN.rounds);

		
		Double error=trainNN(portfolioNN, inputArray, inputData, idealData, idealWeights, trainingRounds);

		portfolioNN.saveNN("portfolioNNwithValdiation-Learning-"+learningRate+"-Rounds-" + portfolioNN.rounds+"Error"+error);

		long endTime = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("total time:" + totalTime);
	}

	public static double trainNN(MLNN portfolioNN, double[][] inputArray, double[][][] inputArrayFull, double[][][] idealData, double[][] idealWeights, int rounds) {
		double sharp=0.0;
		for (int i = 0; i < rounds; i++) {
			sharp=0.0;
			for (int j = 0; j < 150; j++) {
				sharp+=sharpCalc(portfolioNN.forwardProp(inputArray[j], idealWeights[j]), inputArrayFull[j], idealData[j][11][11]);
				portfolioNN.backProp(idealWeights[j]);
			}
			System.out.print("ROUND:" + portfolioNN.totalRounds + "\tTRAINING ERROR:" + portfolioNN.getError() + "\tTRAINING Sharp Average:" + sharp/150);
			sharp=0.0;
			for (int j = 150; j < 250; j++) {
				sharp+=sharpCalc(portfolioNN.forwardProp(inputArray[j], idealWeights[j]), inputArrayFull[j], idealData[j][11][11]);
			}
			System.out.println("\tVALIDATION ERROR:" + portfolioNN.getError()+ "\tVALIDATION Sharp Average:" + sharp/100);
		}
		return sharp/100;
	}

	public static double[][] read0(String inputFile) throws IOException {
		double[][] data = new double[250][10];
		File inputWorkbook = new File(inputFile);
		Workbook w;
		try {
			w = Workbook.getWorkbook(inputWorkbook);
			// Get the first sheet
			Sheet sheet = w.getSheet(0);

			for (int i = 1; i <= 250; i++) {
				for (int j = 1; j <= 10; j++) {
					Cell cell = sheet.getCell(j - 1, i - 1);
					data[i - 1][j - 1] = Double.parseDouble(cell.getContents());
				}
			}
		} catch (BiffException e) {
			e.printStackTrace();
		}
		return data;
	}

	public static double[][] read1(String inputFile) throws IOException {
		double[][] data = null;
		File inputWorkbook = new File(inputFile);
		Workbook w;

		try {
			w = Workbook.getWorkbook(inputWorkbook);

			// Get the first sheet
			Sheet sheet = w.getSheet(0);
			data = new double[60][14];

			for (int j = 1; j < 15; j++) {
				for (int i = 1; i < 61; i++) {
					Cell cell = sheet.getCell(j, i);
					data[i - 1][j - 1] = Double.parseDouble(cell.getContents());
				}
			}

		} catch (BiffException e) {
			e.printStackTrace();
		}
		return data;
	}

	public static double[][] read2(String inputFile) throws IOException {
		double[][] data = null;
		File inputWorkbook = new File(inputFile);
		Workbook w;

		try {
			w = Workbook.getWorkbook(inputWorkbook);

			// Get the first sheet
			Sheet sheet = w.getSheet(1);
			data = new double[11][10];
			for (int j = 1; j < 11; j++) {
				for (int i = 1; i < 12; i++) {
					Cell cell = sheet.getCell(j, i);
					data[i - 1][j - 1] = Double.parseDouble(cell.getContents());
				}
			}

		} catch (BiffException e) {
			e.printStackTrace();
		}
		return data;
	}

	public static double[][] read3(String inputFile) throws IOException {
		double[][] data = null;
		File inputWorkbook = new File(inputFile);
		Workbook w;

		try {
			w = Workbook.getWorkbook(inputWorkbook);

			// Get the first sheet
			Sheet sheet = w.getSheet(2);
			data = new double[12][12];
			for (int j = 1; j < 13; j++) {
				for (int i = 1; i < 13; i++) {
					Cell cell = sheet.getCell(j, i);
					data[i - 1][j - 1] = Double.parseDouble(cell.getContents());
				}
			}

		} catch (BiffException e) {
			e.printStackTrace();
		}
		return data;
	}

	public static double read4(String inputFile) throws IOException {
		File inputWorkbook = new File(inputFile);
		Workbook w;

		try {
			w = Workbook.getWorkbook(inputWorkbook);
			// Get the first sheet
			Sheet sheet = w.getSheet(2);			
			Cell cell = sheet.getCell(12, 12);
			return Double.parseDouble(cell.getContents());	
		} catch (BiffException e) {
			e.printStackTrace();
		}
		return 0.0;
	}
	
	public static void print3d(double[][][] input) {
		for (int i = 0; i < 1; i++) {
			for (int j = 0; j < input[0].length; j++) {

				for (int k = 0; k < input[i][j].length; k++) {

					System.out.print(input[i][j][k] + "\t");
				}
				System.out.println();
			}
		}
	}

	public static double sharpCalc(double[] weights, double[][] idealData, double riskFreeRate) {
		// calulating ideal output array from ideal data

		// mean of portfolios 
		double[] idealMean = new double[10];
		for (int j = 1; j <= 10; j++) {
			double sum = 0;
			for (int k = 1; k <= 12; k++) {
				sum += idealData[k - 1][j];
			}
			idealMean[j - 1] = sum / 12;
		} // std of portfolio
		
		double[] idealSTD = new double[10];
		for (int j = 1; j <= 10; j++) {
			double sum = 0;
			for (int k = 1; k <= 12; k++) {
				sum += (idealData[k - 1][j] - idealMean[j - 1])
						* (idealData[k - 1][j] - idealMean[j - 1]);
			}
			idealSTD[j - 1] = Math.sqrt(sum / 11);
		} // correlation
		double[][] corrIdeal = new double[10][10];

		double[] deltaSquaredSum = new double[10];
		double[][] delta = new double[10][12];

		for (int cola = 1; cola <= 10; cola++) {
			deltaSquaredSum[cola - 1] = 0;
			for (int month = 1; month <= 12; month++) {
				delta[cola - 1][month - 1] = idealData[month - 1][cola] - idealMean[cola - 1];
				deltaSquaredSum[cola - 1] += delta[cola - 1][month - 1] * delta[cola - 1][month - 1];
			}
		}

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				double abSUM = 0;
				for (int month = 1; month <= 12; month++) {
					abSUM += delta[i][month - 1] * delta[j][month - 1];
				}
				corrIdeal[j][i] = abSUM / Math.sqrt(deltaSquaredSum[i] * deltaSquaredSum[j]);
			}
		}
		double[][] stdMatrix = new double[10][10];

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				stdMatrix[i][j] = idealSTD[i] * idealSTD[j];
			}
		}

		// final scalar matrix 

		double[][] CorrSTDMatrix = new double[10][10];

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				CorrSTDMatrix[i][j] = stdMatrix[i][j] * corrIdeal[j][i];
			}
		}

		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				//System.out.print(CorrSTDMatrix[j][i] + "\t");
			}
			//System.out.println();
		}
		
		//weight matrix.. or variance calc
		double var=0.0;
		for (int i = 0; i < 10; i++) {
			for (int j = 0; j < 10; j++) {
				//System.out.print(weights[i] * weights[j]*CorrSTDMatrix[i][j] + "\t");

				var+= weights[i] * weights[j]*CorrSTDMatrix[i][j];
			}
			//System.out.println();

		}		
		
		//System.out.println("var"+ var);
		
		double std=Math.sqrt(var);
		
		//System.out.println("std"+ std);
		
		double totalReturn=0.0;
		for (int i=0; i<10; i++){
			totalReturn +=weights[i]*idealMean[i];
			//System.out.print(idealMean[i] + "\t" + weights[i] + "\t");
		}
		//System.out.println(totalReturn);
		
		return (totalReturn-riskFreeRate)/std;
			
	}
}

/*
 * 
 * // constraints Collection constraints = new ArrayList(); // first constraint,
 *
 * 
 * they all add up to 1 double[] c1 = new double[100];for( int i = 0;i<10;*i++)
 * { c1[10 * i + i] = CorrSTDMatrix[i][i]; }**constraints.add(new
 * LinearConstraint(c1,Relationship.EQ,1.00));** // second constraint, //
 * optimization begins // MultivariateMultiStartOptimizer(MultivariateOptimizer
 * optimizer, int starts, // RandomVectorGenerator generator) /*
 **
 * MultivariateMultiStartOptimizer MMSO = new
 * MultivariateMultiStartOptimizer(MO, 5, new RandomVectorGenerator generator)**
 * MultivariateFunction sharp = new MultivariateFunction() {
 * 
 * @Override public double sharp(double[][] CorrSTDMatrix, double[] weights,
 * double[] mean, double riskfree) { double[][] weightMatrix=new double[10][10];
 * double sum=0; double returnSum=0; for (int i=0; i<10; i++) {
 * returnSum+=mean[i]*weights[i]; for (int j=0; j<10; j++) {
 * sum+=CorrSTDMatrix[i][j]*weights[j]*weights[i]; } } return
 * (returnSum-riskfree)/Math.sqrt(sum); } };
 * 
 */
/*
 * public double sharp(double[][] CorrSTDMatrix, double[] weights, double[]
 * mean, double riskfree) { double[][] weightMatrix=new double[10][10]; double
 * sum=0; double returnSum=0; for (int i=0; i<10; i++) {
 * returnSum+=mean[i]*weights[i]; for (int j=0; j<10; j++) {
 * sum+=CorrSTDMatrix[i][j]*weights[j]*weights[i]; } } return
 * (returnSum-riskfree)/Math.sqrt(sum); }
 * 
 * /* Model model = new Model("SharpMax");
 * 
 * RealVar[] weights = new RealVar[10]; for (int i=0; i<10; i++) {
 * weights[i]=model.realVar(Integer.toString(i), -2.00d, 2.00d, .01); } RealVar
 * sharp=model.realVar("Sharp", 0.00d, 5.00d, .01);
 * 
 * model.realIbexGenericConstraint(
 * "{0}+{1}+{2}+{3}+{4}+{5}+{6}+{7}+{8}+{9}<=1.00",
 * weights[0],weights[1],weights[2],weights[3],weights[4],weights[5],weights[6],
 * weights[7],weights[8],weights[9]).post();
 * 
 * //function is for each pair of funds, multiply the weights, multiply the
 * corresonding correlation matrix values and SD matrix values String
 * sharpFunc=new String("");
 * 
 * for (int i=0; i<10; i++) { for (int j=0; j<10; j++) {
 * sharpFunc.concat("{"+i+"}*"+"{"+j+"}*"+Double.toString(weightMatrix[i][j])+
 * "*"); } }
 * 
 * System.out.println(sharpFunc);
 * 
 * //constraints
 * 
 * model.realIbexGenericConstraint(
 * "{0}+{1}+{2}+{3}+{4}+{5}+{6}+{7}+{8}+{9}={Sharp}",
 * weights[0],weights[1],weights[2],weights[3],weights[4],weights[5],weights[6],
 * weights[7],weights[8],weights[9]).post();
 * 
 * 
 * 
 * 
 * 
 * 
 * // Model objective function 3X + 4Y //IntVar OBJ = model.intVar("objective",
 * 0, 999); //model.scalar(new IntVar[]{X,Y}, new int[]{3,4}, OBJ)).post(); //
 * Specify objective //model.setObjective(Model.MAXIMIZE, OBJ); // Compute
 * optimum model.getSolver().solve(); Solver solver = model.getSolver();
 * Solution solution = solver.findSolution();
 * 
 * } catch (IOException e) { // TODO Auto-generated catch block
 * e.printStackTrace(); } }
 */