package cure;
import java.util.StringTokenizer;

/**
 * Represents a Point Class. Also stores the KD Tree index for search.
 */
public class Point {
	public double x;
	public double y;
	public int index;
	
	public Point() {
	}
	
	public Point(double x, double y, int index) {
		this.x = x;
		this.y = y;
		this.index = index;
	}
	
	public Point(Point point) {
		this.x = point.x;
		this.y = point.y;
	}
	
	public double[] toDouble() {
		double[] xy = {x,y};
		return xy; 
	}
	
	public static Point parseString(String str) {
		Point point = new Point();
		StringTokenizer st = new StringTokenizer(str);
		return point;
	}
	
	/**
	 * Calculates the Euclidean Distance from a Point t
	 */
	public double calcDistanceFromPoint(Point t) {
		return Math.sqrt(Math.pow(x-t.x, 2) + Math.pow(y-t.y, 2));
	}
	
	public String toString() {
		return "{" + x + "," + y + "}";
	}
	
	public boolean equals(Point t) {
		return (x == t.x) && (y == t.y);
	}
}
