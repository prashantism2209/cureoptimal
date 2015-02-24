package cure;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;



public class KDClusterNode {
	static int noOfDimension;
	static int clustercount;
	public int totalNumberOfPoints;
	
	public PointsNDim [] kdpoints;
	public PointsNDim medianCluster;
	public	int level;		
	public	double [] minDimesion; //For density calculation
	public	double [] maxDimesion; //For density Calculation
	public	int pointsInCluster;	
	public KDClusterNode left;
	public KDClusterNode right;
	public KDClusterNode parent;
	public	double density;
	public boolean isLeaf;
	static double [] norMin;
	static double [] norMax;
	static double [] norRange;
	KDClusterNode(int size){
		kdpoints=new PointsNDim[size];
		minDimesion= new double[size];
		maxDimesion= new double[size];
		medianCluster=new PointsNDim(size);
		left=null;
		right=null;
	}
	//calculate values of median for cluster and minDimension and maxDimesion
	/*void calculateMedianMinMax(){
		try{
			double [] dim=new double[noOfDimension];
			double  xmax=this.kdpoints[0].x;
			double  xmin=this.kdpoints[0].y;
			double  ymax=this.kdpoints[0].x;
			double  ymin=this.kdpoints[0].y;
			double [] sum=new double[noOfDimension];
				
			for(int i=0;i<this.pointsInCluster;i++){
				Point point=this.kdpoints[i];
				sum[0]+=point.x;
				sum[1]+=point.y;
				if(xmax<point.x) xmax=point.x;
				if(ymax<point.y) ymax=point.y;
				if(xmin>point.x) xmin=point.x;
				if(ymin>point.y) ymin=point.y;
			}
			//System.out.println(" sum[0]="+sum[0]+" sum[1]"+sum[1]+" points="+this.pointsInCluster+" "+sum[0]/this.pointsInCluster);
			
			if(this.pointsInCluster>0){
				double s1=(double)(sum[0]/this.pointsInCluster);
				double s2=(double)(sum[1]/this.pointsInCluster);
			//	System.out.println("medians "+s1+"  "+s2);
				this.medianCluster.x=s1;
				this.medianCluster.y=s2;
				this.minDimesion[0]=xmin;
				this.minDimesion[1]=ymin;
				this.maxDimesion[0]=xmax;
				this.maxDimesion[1]=ymax;
			//	System.out.println("x "+this.medianCluster.x+ "  y="+this.medianCluster.y+" sum[0]="+sum[0]+" sum[1]"+sum[1]+" points="+this.pointsInCluster);
				
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}*/
	void calculateMedianMinMax(){
		try{
			double [] dim=new double[noOfDimension];
			double [] dmax=new double[noOfDimension];
			double [] dmin=new double[noOfDimension];			
			double [] sum=new double[noOfDimension];
				
			for(int j=0;j<noOfDimension;j++){
				dmax[j]=-1000000.0;
				dmin[j]=10000000.0;
				sum[j]=0.0;				
			}
			
			for(int i=0;i<this.pointsInCluster;i++){
				
				for(int j=0;j<noOfDimension;j++){
					sum[j]+=this.kdpoints[i].coord[j];
					if(dmax[j]<this.kdpoints[i].coord[j]){
						dmax[j]=this.kdpoints[i].coord[j];
					}
					if(dmin[j]>this.kdpoints[i].coord[j]){
						dmin[j]=this.kdpoints[i].coord[j];
					}
				}
				
			}
			//System.out.println(" sum[0]="+sum[0]+" sum[1]"+sum[1]+" points="+this.pointsInCluster+" "+sum[0]/this.pointsInCluster);
			
			if(this.pointsInCluster>0){
				
				for(int j=0;j<noOfDimension;j++){
					
					this.medianCluster.coord[j]=sum[j]/this.pointsInCluster;
					this.minDimesion[j]=dmin[j];
					this.maxDimesion[j]=dmax[j];
				}
			//	System.out.println("medians "+s1+"  "+s2);
			
			//	System.out.println("x "+this.medianCluster.x+ "  y="+this.medianCluster.y+" sum[0]="+sum[0]+" sum[1]"+sum[1]+" points="+this.pointsInCluster);
				
			}
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	void calculateDensity(){
		double volume=1;
		for(int j=0;j<noOfDimension;j++){
			volume=volume*(this.maxDimesion[j]-this.minDimesion[j]);
			//volume=volume*(this.maxDimesion[j]-this.minDimesion[j])/(this.maxDimesion[j]*this.minDimesion[j]);
		}
		this.density=this.pointsInCluster/volume;
	}
	
	static void buildTree(KDClusterNode root, int level){
		
		KDClusterNode left=new KDClusterNode(root.pointsInCluster);
		KDClusterNode right=new KDClusterNode(root.pointsInCluster);
		//setting parent node
		left.parent=root;
		right.parent=root;
		root.left=left;
		root.right=right;
		
		int leftClusterPoint=0;
		int rightClusterPoint=0;	
		
		root.removeNoise(root);
		
		for(int i=0;i<root.pointsInCluster;i++){
			PointsNDim point=root.kdpoints[i];
			
				if(point.coord[level]<root.medianCluster.coord[level]){					
					left.kdpoints[leftClusterPoint]=point;
					leftClusterPoint++;		
				}else{
					right.kdpoints[rightClusterPoint]=point;
					rightClusterPoint++;					
				}
		
		}
		level=(level+1)%noOfDimension;
		
		left.pointsInCluster=leftClusterPoint;
		right.pointsInCluster=rightClusterPoint;
		
		left.calculateMedianMinMax();
		right.calculateMedianMinMax();
		
		left.calculateDensity();
		right.calculateDensity();
		
		left.isLeaf=true;
		int treeflag=0;
		if(left.pointsInCluster>5&&left.density>root.density){
			left.isLeaf=false;
			treeflag=1;
				KDClusterNode.buildTree(left,level);
			
		}else{
			left.isLeaf=true;
		}
		
		right.isLeaf=true;
		if(right.pointsInCluster>5&&right.density>right.density){
			right.isLeaf=false;
			treeflag=1;
				KDClusterNode.buildTree(right,level);
			
		}else{
			right.isLeaf=true;
		}
		if(treeflag==0){
			root.isLeaf=true;
			root.left=null;
			root.right=null;
		}
	}
	
  public void removeNoise(KDClusterNode root){
	  	double sum=0.0;
		for(int j=0;j<noOfDimension;j++){
			sum=sum+((this.maxDimesion[j]-this.minDimesion[j])*(this.maxDimesion[j]-this.minDimesion[j]));
		}
		  double thresholdDistance= 1*Math.sqrt(sum);
		  System.out.println("Threshold Distance="+thresholdDistance);
		  System.out.println("Cluster");
		  for(int i=0;i<root.pointsInCluster;i++){
				PointsNDim point=root.kdpoints[i];
				
					if(PointsNDim.sqrdist(point,root.medianCluster)>thresholdDistance){
						//System.out.println("Out lier"+point.toString());
					}
				
			}
				  
				 
	  }
	
	 public static void inOrder (KDClusterNode root)
	 {
		  if(root == null) return;
		 if(root.isLeaf==true)
		 { 
			 clustercount++;	
			 root.logClusterParallel(root, "iris_1.csv","a"+clustercount);
		 }
		  inOrder( root.left );		 	  		  
		  inOrder( root.right ); 		  
	}
	  public static void merge(KDClusterNode root){
		  if(root == null) return;
		  
		  if(root.isLeaf==true)
		  {
			  System.out.println(root.pointsInCluster);	
			  
		  }
		  else{
			  merge( root.left );	
			  merge( root.right ); 
		  }
		 
	  }
	  
	  private static void closeWriterHandle(BufferedWriter out) {
			try {
				out.flush();
				out.close();
			} catch(Exception e) {
				//debug(e);
			}
		}
		
	  private static BufferedWriter getWriterHandle(String filename) {
			BufferedWriter out = null;
			try {
				FileWriter fw = new FileWriter(filename, true);
				out = new BufferedWriter(fw);
			} catch(Exception e) {
				//debug(e);
			}
			return out;
		}
		private void logCluster(KDClusterNode node, String filename) {
			BufferedWriter out = getWriterHandle(filename);
			try {
				out.write("name,a,b,c,d" + "\n");
				for(int j=0; j<node.pointsInCluster; j++) {
					PointsNDim p = (PointsNDim)node.kdpoints[j];
					out.write(""+p.index);
					for(int t=0;t<KDClusterNode.noOfDimension-1;t++){
						out.write("," + p.coord[t]);
					}
					
					out.write("," + p.coord[KDClusterNode.noOfDimension-1] + "\n");
					
				}
				
			} catch(Exception e){
				//debug(e);
			}
			closeWriterHandle(out);
		}
		private void logClusterParallel(KDClusterNode node, String filename,String name) {
			BufferedWriter out = getWriterHandle(filename);
			try {
				//out.write("name,a,b,c,d" + "\n");
				for(int j=0; j<node.pointsInCluster; j++) {
					PointsNDim p = (PointsNDim)node.kdpoints[j];
					out.write(""+name);
					for(int t=0;t<KDClusterNode.noOfDimension-1;t++){
						out.write("," + p.coord[t]);
					}
					
					out.write("," + p.coord[KDClusterNode.noOfDimension-1] + "\n");
					
				}
				
			} catch(Exception e){
				//debug(e);
			}
			closeWriterHandle(out);
		}
		private void logInitializePallel(String filename){
			BufferedWriter out = getWriterHandle(filename);
			try {
				out.write("cluster,a,b,c,d" + "\n");
				//out.write("cluster,a,b,c" + "\n");
								
			} catch(Exception e){
				//debug(e);
			}
			closeWriterHandle(out);
		}
	  protected static void readDataPoints(KDClusterNode KDNode) {
			int pointIndex = 0;
			FileReader fr = null;
			System.out.println("readdata");
			//File fl=new File("spaeth2_05.txt");
			File fl=new File("iris.txt");
			//File fl=new File("3dimension.txt");
			//File fl=new File("set.txt");
			try {
				//fr = new FileReader(dataFile);
				fr = new FileReader(fl); 
				BufferedReader in = new BufferedReader(fr);
				String data = in.readLine();
				//System.out.println("print first data="+data);
				int cnt=0;
				while(data != null) {
					StringTokenizer st = new StringTokenizer(data);
					for(int j=0;j<KDClusterNode.noOfDimension;j++){
						double x = Double.parseDouble(st.nextToken());
						KDNode.kdpoints[cnt].coord[j]=x;
						KDNode.kdpoints[cnt].index=cnt;
					}
					
					System.out.println(KDNode.kdpoints[cnt].toString());
					cnt++;
					
					data = in.readLine();
				}
				in.close();
			} catch(Exception e){
				e.getStackTrace();
				}
			}
	  protected static void normalizePoints(KDClusterNode node)
	  {
		  KDClusterNode.norMax= new double[KDClusterNode.noOfDimension];
		  KDClusterNode.norMin= new double[KDClusterNode.noOfDimension];
		  KDClusterNode.norRange= new double[KDClusterNode.noOfDimension];
			for(int t=0;t<KDClusterNode.noOfDimension;t++)
			{
				
				
				KDClusterNode.norMax[t]=-100000000;
				KDClusterNode.norMin[t]=100000000;
				
			}
		  for(int j=0; j<node.pointsInCluster; j++) {
				PointsNDim p = (PointsNDim)node.kdpoints[j];
				
				for(int t=0;t<KDClusterNode.noOfDimension;t++){
					
					if(p.coord[t]>KDClusterNode.norMax[t])
						KDClusterNode.norMax[t]=p.coord[t];
					if(p.coord[t]<KDClusterNode.norMin[t])
						KDClusterNode.norMin[t]=p.coord[t];
					KDClusterNode.norRange[t]=KDClusterNode.norMax[t]-KDClusterNode.norMin[t];
				}
		  }
				for(int j=0; j<node.pointsInCluster; j++)
				{
					PointsNDim p = (PointsNDim)node.kdpoints[j];
					System.out.println("");
					for(int t=0;t<KDClusterNode.noOfDimension;t++)
					{
						
						p.coord[t]=(p.coord[t]-KDClusterNode.norMin[t])/KDClusterNode.norRange[t];
						System.out.print(" "+p.coord[t]);
					}
				
				} 
		  
	  }
	  protected static void denormalizePoints(KDClusterNode node)
	  {
		 
		  
				for(int j=0; j<node.pointsInCluster; j++)
				{
					PointsNDim p = (PointsNDim)node.kdpoints[j];
					System.out.println("");
					for(int t=0;t<KDClusterNode.noOfDimension;t++)
					{
						p.coord[t]=(p.coord[t]*KDClusterNode.norRange[t])+KDClusterNode.norMin[t];
						System.out.print(" "+p.coord[t]);
					}
				
				} 
		  
	  }
	  //Main function to be removed later
	  public static void main(String[] args){
		   KDClusterNode KDNode=new KDClusterNode(150);
		  	KDNode.pointsInCluster=150;
			KDClusterNode.noOfDimension=4;
			
			for(int j=0;j<KDNode.pointsInCluster;j++){
				KDNode.kdpoints[j]=new PointsNDim(KDClusterNode.noOfDimension);
			}
		//	KDNode.kdpoints=dataPoints;
			KDClusterNode.readDataPoints(KDNode);
			//normalize
			KDClusterNode.normalizePoints(KDNode);
			KDNode.parent=null;
			KDNode.left=null;
			KDNode.right=null;
			
			System.out.println("pointss=="+KDNode.pointsInCluster);
			KDNode.calculateMedianMinMax();
			KDNode.calculateDensity();
			if(KDNode.pointsInCluster>10){
				KDNode.isLeaf=false;
				KDClusterNode.buildTree(KDNode,0);
			}else{
				KDNode.isLeaf=true;
			}
			//denormalize
			System.out.println("Denormalize");
			KDClusterNode.denormalizePoints(KDNode);
			KDClusterNode.clustercount=0;
			KDNode.logInitializePallel("iris_1.csv");
			KDClusterNode.inOrder(KDNode);
			KDNode.logCluster(KDNode, "cars.csv");
			System.out.println("no of Cluster=="+KDClusterNode.clustercount);
		}
	 
}
