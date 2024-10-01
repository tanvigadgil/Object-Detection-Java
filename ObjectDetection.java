import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.stream.Collectors;
import java.util.List;
import java.util.Map;
import javax.swing.*;

public class ObjectDetection {
    private static final int HUE_RANGE = 360;
    private static final double NORMALIZED_HUE_THRESHOLD = 3.0 / HUE_RANGE;
    private static final int GRID_SIZE = 15;
	public static final Color COLOUR = new Color(102, 0, 153); // Purple

    JFrame frame;
    JLabel lbIm1;
    BufferedImage imgIpt, imgObj;
    int width = 640;
    int height = 480;
    List<int[]> objectPixelCoordinates = new ArrayList<>();

	// ********************** Point class **********************
	private static class Point {
        int x;
        int y;

        Point(int x, int y) {
            this.x = x;
            this.y = y;
        }

        Point(int[] coordinates) {
            this.x = coordinates[0];
            this.y = coordinates[1];
        }
    }

	// ********************** GridCell class **********************
	private static class GridCell {
        private final int x;
        private final int y;

        public GridCell(int x, int y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o)
                return true;
            if (o == null || getClass() != o.getClass())
                return false;
            GridCell gridCell = (GridCell) o;
            return x == gridCell.x && y == gridCell.y;
        }

        @Override
        public int hashCode() {
            return 31 * x + y;
        }
    }

	// ********************** Functions Begin **********************
	// Euclidean distance between two points
	private static double distance(Point p1, Point p2) {
        int dx = p1.x - p2.x;
        int dy = p1.y - p2.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

	// RGB to HSV conversion
    private double[] rgbToHsv(double r, double g, double b) {
        r = r / 255.0;
        g = g / 255.0;
        b = b / 255.0;

        double cmax = Math.max(r, Math.max(g, b));
        double cmin = Math.min(r, Math.min(g, b));

        double diff = cmax - cmin;
        double h = -1, s = -1;

        if (cmax == cmin) {
            h = 0;
		} else if (cmax == r) {
            h = (60 * ((g - b) / diff) + 360) % 360;
		} else if (cmax == g) {
            h = (60 * ((b - r) / diff) + 120) % 360;
		} else if (cmax == b) {
            h = (60 * ((r - g) / diff) + 240) % 360;
		}

        s = (cmax == 0) ? 0 : (diff / cmax) * 100;

        double v = cmax * 100;
        return new double[] { h, s, v };
    }

	// Read Image (Original)
    private void readImageRGB(int width, int height, String imgPath, BufferedImage imgHSV) {
        try {
            int frameLength = width * height * 3;
            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            raf.seek(0);
            long len = frameLength;
            byte[] bytes = new byte[(int) len];

            raf.read(bytes);

            int ind = 0;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    byte a = 0;
                    int r = bytes[ind] & 0xFF;
                    int g = bytes[ind + height * width] & 0xFF;
                    int b = bytes[ind + height * width * 2] & 0xFF;

                    double[] hsv = rgbToHsv(r, g, b);

                    int hsvColor = Color.HSBtoRGB((float) (hsv[0] / 360.0), (float) (hsv[1] / 100.0),
                            (float) (hsv[2] / 100.0));
                    imgHSV.setRGB(x, y, hsvColor);
                    ind++;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

	// Compute histogram and store coordinates of input image
    private List<int[]> inputImageHist(BufferedImage imgHSV) {

		List<int[]> inputCoordinates = new ArrayList<>();
        int[] hist = new int[360];
		int maxValue;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = imgHSV.getRGB(x, y);
                float[] hsv = Color.RGBtoHSB((rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF, null);
                int hue = (int) (hsv[0] * 360);
                hist[hue]++;
                inputCoordinates.add(new int[] { x, y, hue });
            }
        }

        maxValue = Arrays.stream(hist).max().orElse(1);

        for (int i = 0; i < hist.length; i++) {
            hist[i] = (int) (hist[i] * (double) height / maxValue);
        }
        return inputCoordinates;
    }

	// Compute histogram and store coordinates of object image
	// Remove Green chroma background
    private List<int[]> objectHist(BufferedImage imgHSV) {

        List<int[]> objectCoordinates = new ArrayList<>();
        int[] hist = new int[360];
		int maxValue;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = imgHSV.getRGB(x, y);
                float[] hsv = Color.RGBtoHSB((rgb >> 16) & 0xFF, (rgb >> 8) & 0xFF, rgb & 0xFF, null);
                int hue = (int) (hsv[0] * 360);

                if (hue != 120) {
                    hist[hue]++;
                    objectCoordinates.add(new int[] { x, y, hue });
                }
            }
        }

        maxValue = Arrays.stream(hist).max().orElse(1);
        for (int i = 0; i < hist.length; i++) {
            hist[i] = (int) (hist[i] * (double) height / maxValue);
        }

        return objectCoordinates;
    }

	private static Map<GridCell, List<int[]>> buildGrid(List<int[]> coordinates) {
        Map<GridCell, List<int[]>> grid = new HashMap<>();
        for (int[] coordinate : coordinates) {
            GridCell cell = getGridCell(coordinate);
            grid.computeIfAbsent(cell, k -> new ArrayList<>()).add(coordinate);
        }
        return grid;
    }

	private static GridCell getGridCell(int[] coordinate) {
        int x = coordinate[0] / GRID_SIZE;
        int y = coordinate[1] / GRID_SIZE;
        return new GridCell(x, y);
    }

	// Calculate matching hues of coordinates between input and object images
    private static List<int[]> calcMatchingCoords(List<int[]> inputCoordinates, List<int[]> objectCoordinates) {

        List<int[]> matchingCoordinates = new ArrayList<>();
        Map<GridCell, List<int[]>> objectGrid = buildGrid(objectCoordinates);

        for (int[] inputCoordinate : inputCoordinates) {
            GridCell cell = getGridCell(inputCoordinate);
            boolean matchFound = false;

            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    GridCell neighborCell = new GridCell(cell.x + dx, cell.y + dy);
                    List<int[]> matches = objectGrid.get(neighborCell);

                    if (matches != null) {
                        for (int[] match : matches) {
                            double normInputHue = inputCoordinate[2] / (double) HUE_RANGE;
                            double normObjectHue = match[2] / (double) HUE_RANGE;
                            double hueDiff = Math.abs(normObjectHue - normInputHue);

                            double spatialDiff = Math.sqrt(Math.pow(match[0] - inputCoordinate[0], 2) + Math.pow(match[1] - inputCoordinate[1], 2));
                            double totalDiff = 0.7 * hueDiff + 0.9 * spatialDiff;

                            if (totalDiff <= NORMALIZED_HUE_THRESHOLD) {
                                matchFound = true;
                                break;
                            }
                        }
                    }
                    if (matchFound) {
                        break;
                    }
                }

                if (matchFound) {
                    break;
                }
            }
            if (matchFound) {
                matchingCoordinates.add(inputCoordinate);
            }
        }

        return matchingCoordinates;
    }


	// Check and merge overlapping bounding boxes
	private List<Rectangle> mergeOverlappingBoundingBoxes(List<Rectangle> boundingBoxes) {
		List<Rectangle> mergedBoundingBoxes = new ArrayList<>();
	
		for (Rectangle currentBoundingBox : boundingBoxes) {
			boolean merged = false;
			for (int i = 0; i < mergedBoundingBoxes.size(); i++) {
				Rectangle mergedBoundingBox = mergedBoundingBoxes.get(i);
				
				Rectangle intersection = currentBoundingBox.intersection(mergedBoundingBox);
				
				if (!intersection.isEmpty()) {
					mergedBoundingBox.add(currentBoundingBox);
					merged = true;
					break;
				}
			}
			
			if (!merged) {
				mergedBoundingBoxes.add(currentBoundingBox);
			}
		}
	
		return mergedBoundingBoxes;
	}
	
	// Draw bounding boxes around clusters and display the name of the object
    private void drawBoundingBox(BufferedImage img, List<List<Point>> clusters, int borderThickness, String objectName) {
        if (clusters.isEmpty()) {
            System.out.println("Object not found.");
            return;
        }

		List<Rectangle> boundingBoxes = new ArrayList<>();
		List<Rectangle> mergedBoundingBoxes;

        Graphics2D g = img.createGraphics();
        g.setColor(COLOUR);

		Font font = new Font("Arial", Font.BOLD, 15);
    	g.setFont(font);

        for (List<Point> cluster : clusters) {
            int minX = cluster.stream().mapToInt(point -> point.x).min().orElse(0);
            int minY = cluster.stream().mapToInt(point -> point.y).min().orElse(0);
            int maxX = cluster.stream().mapToInt(point -> point.x).max().orElse(0);
            int maxY = cluster.stream().mapToInt(point -> point.y).max().orElse(0);

            minX = Math.max(0, minX - borderThickness);
            minY = Math.max(0, minY - borderThickness);
            maxX = Math.min(img.getWidth() - 1, maxX + borderThickness);
            maxY = Math.min(img.getHeight() - 1, maxY + borderThickness);

            int width = maxX - minX + 1;
            int height = maxY - minY + 1;

            for (int i = 0; i < borderThickness; i++) {
				boundingBoxes.add(new Rectangle(minX + i, minY + i, width - 2 * i, height - 2 * i));
            }
        }

		mergedBoundingBoxes= mergeOverlappingBoundingBoxes(boundingBoxes);

		for (int i = 0; i < borderThickness; i++) {
			for (Rectangle mergedBoundingBox : mergedBoundingBoxes) {
				g.drawRect(mergedBoundingBox.x + i, mergedBoundingBox.y + i, mergedBoundingBox.width - 2 * i, mergedBoundingBox.height - 2 * i);

				g.drawString(objectName, mergedBoundingBox.x + borderThickness, mergedBoundingBox.y + mergedBoundingBox.height + 12);
			}
		}
		
        g.dispose();
    }

	// Filter out small objects from the Input image
    private static List<int[]> filterSmallObjects(List<int[]> matchingCoordinates, int minSize, double maxDistance) {
        List<Point> points = matchingCoordinates.stream().map(coord -> new Point(coord[0], coord[1])).collect(Collectors.toList());

        List<List<Point>> clusters = clusterPoints(points, maxDistance);
        List<int[]> filteredCoordinates = clusters.stream().filter(cluster -> cluster.size() >= minSize).flatMap(Collection::stream).map(point -> new int[] { point.x, point.y }).collect(Collectors.toList());

        return filteredCoordinates;
    }

	// Compute the clusters
    private static List<List<Point>> clusterPoints(List<Point> points, double maxDistance) {
        List<List<Point>> clusters = new ArrayList<>();

        for (Point point : points) {
            boolean assigned = false;

            for (List<Point> cluster : clusters) {
                Point clusterPoint = cluster.get(0);

                if (distance(point, clusterPoint) <= maxDistance) {
                    cluster.add(point);
                    assigned = true;
                    break;
                }
            }

            if (!assigned) {
                List<Point> newCluster = new ArrayList<>();
                newCluster.add(point);
                clusters.add(newCluster);
            }
        }

        List<List<Point>> mergedClusters = new ArrayList<>();

        for (List<Point> cluster : clusters) {
            boolean merged = false;

            for (List<Point> mergedCluster : mergedClusters) {
                Point clusterPoint = cluster.get(0);
                Point mergedClusterPoint = mergedCluster.get(0);

                if (distance(clusterPoint, mergedClusterPoint) <= maxDistance) {
                    mergedCluster.addAll(cluster);
                    merged = true;
                    break;
                }
            }

            if (!merged) {
                mergedClusters.add(cluster);
            }
        }

        return mergedClusters;
    }

	// Display the image (Original)
    public void showIms(String[] args) {
        if (args.length < 2) {
            System.out.println("Too few arguments.");
            return;
        }

        imgIpt = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, args[0], imgIpt);

        frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        frame.getContentPane().setLayout(gLayout);

        List<int[]> inputCoords = inputImageHist(imgIpt);

        lbIm1 = new JLabel();
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 1;
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 2;
        frame.getContentPane().add(lbIm1, c);

        frame.pack();
        frame.setVisible(true);
        frame.setSize(650, 480);

        for (int i = 1; i < args.length; i++) {
			String objectName = "Object " + i;
			System.out.println("Processing " + objectName + "...");

            imgObj = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            readImageRGB(width, height, args[i], imgObj);
            List<int[]> objectCoords = objectHist(imgObj);

            List<int[]> matchingCoords = calcMatchingCoords(inputCoords, objectCoords);

            List<int[]> filteredCoords = filterSmallObjects(matchingCoords, 40, 7);

            List<Point> filteredPts = filteredCoords.stream().map(coord -> new Point(coord)).collect(Collectors.toList());

            List<List<Point>> allClusters = clusterPoints(filteredPts, 140);

            drawBoundingBox(imgIpt, allClusters, 3, objectName);
			System.out.println("Done!");

        }

		System.out.println("Object Detection Complete!");
		System.out.println("Displaying image...");
        lbIm1.setIcon(new ImageIcon(imgIpt));
    }

    public static void main(String[] args) {
        ObjectDetection ren = new ObjectDetection();
        ren.showIms(args);
    }
}