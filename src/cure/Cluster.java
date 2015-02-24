package cure;
import java.util.ArrayList;

/**
 * Class Cluster represents a collection of points and its set of representative points.
 * It also stores the distance from its closest neighboring cluster.
 */
public class Cluster {

	public ArrayList rep = new ArrayList();
	public ArrayList pointsInCluster = new ArrayList();
	public double distanceFromClosest = 0;
	public Cluster closestCluster;
	public ArrayList closestClusterRep = new ArrayList();
	
	public double computeDistanceFromCluster(Cluster cluster) {
		double minDistance = 1000000;
		for(int i = 0; i<rep.size(); i++) {
			for(int j = 0; j<cluster.rep.size() ; j++) {
				Point p1 = (Point)rep.get(i);
				Point p2 = (Point)cluster.rep.get(j);
				double distance = p1.calcDistanceFromPoint(p2);
				if(minDistance > distance) minDistance = distance;
			}
		}
		return minDistance;
	}
	
	public double density(Cluster cluster) {
		double minDistance = 1000000;
		Point pi = null;
		Point pj=null;
		Point uij=new Point();
		float fx=0;
		double density=0.0;
		double stdDeviation=0.0;
		for(int i = 0; i<rep.size(); i++) {
			for(int j = 0; j<cluster.rep.size() ; j++) {
				Point p1 = (Point)rep.get(i);
				Point p2 = (Point)cluster.rep.get(j);
				double distance = p1.calcDistanceFromPoint(p2);
				if(minDistance > distance){ 
					minDistance = distance;
					pi=p1;
					pj=p2;
				}
			}
		}
		
		uij.x=(pi.x+pj.x)/2;
		uij.y=(pi.y+pj.y)/2;
	
		stdDeviation=(stdDev()+cluster.stdDev())/2;
		
	//	System.out.println("Inside Density Function stdDeviation="+stdDeviation);
		
		for(int i = 0; i<getClusterSize(); i++) {
			Point p1 = (Point)pointsInCluster.get(i);
			double distance=p1.calcDistanceFromPoint(uij);
		//	System.out.println("Inside Density Function distance="+distance);
			if(distance<=stdDeviation){
				fx=fx+1;
			}
		}
		
		for(int i = 0; i<cluster.getClusterSize(); i++) {
			Point p1 = (Point)cluster.pointsInCluster.get(i);
			double distance=p1.calcDistanceFromPoint(uij);
			//System.out.println("Inside Density Function distance="+distance);
			if(distance<=stdDeviation){
				fx=fx+1;
			}
		}
		density=fx/(getClusterSize()+cluster.getClusterSize());
		
		//System.out.println("Inside Density Function Density="+density);
		return density;
	}
	
	public double intraClusterDistance() {
		double minDistance = 1000000;
		double intraDistance=0.0;
		for(int i = 0; i<getClusterSize(); i++) {
			minDistance = 1000000;
			for(int j = 0; j<rep.size() ; j++) {
				Point p1 = (Point)pointsInCluster.get(i);
				Point p2 = (Point)rep.get(j);
				double distance = p1.calcDistanceFromPoint(p2);
				if(minDistance > distance) minDistance = distance;
			}
			intraDistance+=minDistance;
		}
		
		intraDistance=intraDistance/getClusterSize();
		return intraDistance;
	}
	
	public double intraClusterDensity(double stdDeviation ) {
		
		double intraDensity=0.0;
		//double stdDeviation=stdDev();
		float fx=0;
		
			for(int j = 0; j<rep.size() ; j++) {
				
				for(int i = 0; i<getClusterSize(); i++) {
					
					Point p1 = (Point)pointsInCluster.get(i);
					Point p2 = (Point)rep.get(j);
					double distance = p1.calcDistanceFromPoint(p2);
					if(distance<=stdDeviation) {
						fx=fx+1;
					}
			}
			
		}
			intraDensity=fx/(stdDeviation*rep.size());
			//System.out.println("Intra cluster density="+intraDensity+ "ponts in cluster="+getClusterSize());
			
		return intraDensity;
	}
	
	public double stdDev() {
		double minDistance = 1000000;
		double intraDistance=0.0;
		double sumx=0.0;
		double sumy=0.0;
		double stddevx=0.0;
		double stddevy=0.0;
		double stddev=0.0;
		Point mean=new Point();
		for(int j = 0; j<rep.size() ; j++) {
			
			Point p = (Point)rep.get(j);
			sumx+=p.x;
			sumy+=p.y;
		}
		mean.x=sumx/rep.size();
		mean.y=sumy/rep.size();
		
		double r=0.0;
		for(int j = 0; j<rep.size() ; j++) {
			Point p =  (Point)rep.get(j);
			stddevx=(p.x-mean.x)*(p.x-mean.x);
			stddevy=(p.y-mean.y)*(p.y-mean.y);
			r+=Math.sqrt(stddevx+stddevy);
		}
		double meanr=r/rep.size();
		r=0.0;
		double sdr=0.0;
	for(int j = 0; j<rep.size() ; j++) {
			Point p =  (Point)rep.get(j);
			stddevx+=(p.x-mean.x)*(p.x-mean.x);
			stddevy+=(p.y-mean.y)*(p.y-mean.y);
			r=Math.sqrt(stddevx+stddevy);
			sdr+=(r-meanr)*(r-meanr);
		}
	
	sdr=sdr/rep.size();
	sdr=Math.sqrt(sdr);
		//stddevx/=rep.size();
		//stddevy/=rep.size();
		//stddevx=Math.sqrt(stddevx);
		//stddevx=Math.sqrt(stddevx);
		//stddev=stddevx+stddevy;
		//stddev=Math.sqrt(stddev);
		return sdr;
	}
	
	public int getClusterSize() {
		return pointsInCluster.size();
	}
	
	public ArrayList getPointsInCluster() {
		return pointsInCluster;
	}
}
