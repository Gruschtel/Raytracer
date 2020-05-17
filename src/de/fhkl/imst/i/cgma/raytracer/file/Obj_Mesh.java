package de.fhkl.imst.i.cgma.raytracer.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

public class Obj_Mesh extends RT_Object {

	// ===========================================================
	// Fields
	// ===========================================================

	public int verticesCounter = 0;
	public int textureCounter = 0;
	public int materialCounter = 0;
	public int normalCounter = 0;
	public int facesCounter = 0;
	public int verticesMatCounter = 0;

	public List<String> pMtllib;
	public List<String> pPrimitvType;
	public List<String> pUsemtl;

	public List<float[]> Ka, Kd, Ks; // Ka = ambient, Kd = diffuse, Ks = specular
	public List<float[]> Ns;// ranges between 0 and 1000

	public List<float[]> verticesList; // vertices
	public List<float[]> vertexNormalsList; // vertexNormals
	public List<int[]> facesV; // Faces = triangles
	public List<int[]> facesVN; // Faces = triangles

	public List<int[]> verticesMatList;
	//
	public float[][] materials;
	public int[] materialsN;

	public float[][] vertices;
	public int[] verticesMat;

	public char fgp = 'f'; // flat, gouraud, phong

	// calculated information
	public float[][] vertexNormals; //
	public float[][] vertexColors;

	public int[][] triangles;
	public float[][] triangleNormals;
	public float[][] triangleColors;
	public float[] triangleAreas;

	// ===========================================================
	// Method from SuperClass
	// ===========================================================

	@Override
	public void calcBoundingBox() {
		min[0] = vertices[0][0];
		min[1] = vertices[0][1];
		min[2] = vertices[0][2];

		max[0] = vertices[0][0];
		max[1] = vertices[0][1];
		max[2] = vertices[0][2];

		for (float[] f : vertices) {
			min[0] = min[0] > f[0] ? f[0] : min[0];
			min[1] = min[1] > f[1] ? f[1] : min[1];
			min[2] = min[2] > f[2] ? f[2] : min[2];

			max[0] = max[0] < f[0] ? f[0] : max[0];
			max[1] = max[1] < f[1] ? f[1] : max[1];
			max[2] = max[2] < f[2] ? f[2] : max[2];
		}

		System.out.println(min[0] + " " + min[1] + " " + min[2] + " " + max[0] + " " + max[1] + " " + max[2]);
	}

	@Override
	public String getHeader() {
		return "Obj";
	}

	@Override
	protected void readContent(LineNumberReader br) throws IOException {
		// --------------------------
		// get Values from File
		// --------------------------

		//
		// init Fields
		//
		verticesList = new ArrayList<float[]>();
		verticesMatList = new ArrayList<int[]>();
		vertexNormalsList = new ArrayList<float[]>();
		facesV = new ArrayList<int[]>();
		facesVN = new ArrayList<int[]>();
		Ka = new ArrayList<float[]>();
		Kd = new ArrayList<float[]>();
		Ks = new ArrayList<float[]>();
		Ns = new ArrayList<float[]>();
		pMtllib = new ArrayList<String>();
		pPrimitvType = new ArrayList<String>();
		pUsemtl = new ArrayList<String>();

		//
		// Datei Informationen lesen
		//
		String line;
		while ((line = br.readLine()) != null) {
			line = line.trim();

			StringTokenizer st = new StringTokenizer(line);
			String identifier = st.nextToken().toLowerCase();

			//
			// v: Vertex coordinates
			//
			if (identifier.equals("v")) {
				verticesList.add(new float[] { Float.parseFloat(st.nextToken()), Float.parseFloat(st.nextToken()),
						Float.parseFloat(st.nextToken()) });
				verticesCounter++;
			}

			//
			// vn: Vertex normal
			//
			else if (identifier.equals("vn")) {
				vertexNormalsList.add(new float[] { Float.parseFloat(st.nextToken()), Float.parseFloat(st.nextToken()),
						Float.parseFloat(st.nextToken()) });
				normalCounter++;
			}

			//
			// f: A face definition
			//
			else if (identifier.equals("f")) {
				int tokens = st.countTokens();
				facesV.add(new int[tokens]);
				for (int i = 0; i < tokens; i++) {
					String checkString = st.nextToken();
					facesV.get(facesCounter)[i] = Integer.parseInt(checkString.split("/")[0]);

					// facesVN.get(facesCounter)[i] =
					// Integer.parseInt(st.nextToken().split("//")[1]);
				}
				facesCounter++;
			}

			//
			// mtllib: Name of the MTL file
			//
			else if (identifier.equals("mtllib")) {
				int token = st.countTokens();
				for (int i = 0; i < token; i++) {
					pMtllib.add(st.nextToken());
					materialCounter++;
				}
			}

			//
			// usemtl: Material groups
			//
			else if (identifier.equals("usemtl")) {
				pUsemtl.add(st.nextToken());
				boolean foundMat = false;

				// n
				Ns.add(new float[] { 10f });
				// ambient
				Ka.add(new float[] { 0.15f, 0.15f, 0.15f });
				// diffuse
				Kd.add(new float[] { 0.8f, 0f, 0f });
				// specular
				Ks.add(new float[] { 0.8f, 0.8f, 0.8f });

				for (String fileName : pMtllib) {
					if (new File("data/" + fileName).exists()) {
						for (String lines : readFile("data/" + fileName)) {
							// Nicht benötigte Lines werden gelöscht
							String materialLine = lines.trim();
							if (materialLine.startsWith("#") || materialLine.trim().isEmpty())
								continue;

							StringTokenizer stMaterial = new StringTokenizer(materialLine);
							String identifierMaterial = stMaterial.nextToken().toLowerCase();

							if (identifierMaterial.equals("newmtl")) {
								if (stMaterial.nextToken().equals(pUsemtl.get(textureCounter))) {
									foundMat = true;
								} else {
									foundMat = false;
								}
							}

							if (!foundMat)
								continue;

							//
							// NS: specular exponent
							//
							if (identifierMaterial.equals("ns"))
								Ns.set(textureCounter, new float[] { Float.parseFloat(stMaterial.nextToken()) });
							//
							// Ka: Ambiente [0] = r; [1] = g; [2] = b;
							//
							else if (identifierMaterial.equals("ka")) {

								Ka.set(textureCounter,
										new float[] { Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()) });
								//
								// Kd: Diffuse [0] = r; [1] = g; [2] = b;
								//
							} else if (identifierMaterial.equals("kd")) {
								Kd.set(textureCounter,
										new float[] { Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()) });
								//
								// Ks: Specular [0] = r; [1] = g; [2] = b;
								//
							} else if (identifierMaterial.equals("ks")) {
								Ks.set(textureCounter,
										new float[] { Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()),
												Float.parseFloat(stMaterial.nextToken()) });
							}
						}
					}
				}

				for (int i = verticesMatCounter; i < verticesList.size(); i++) {
					verticesMatList.add(new int[] { textureCounter });
					verticesMatCounter++;
				}
				textureCounter++;
			}

			//
			// Object
			//
			else if (identifier.equals("o")) {
				pPrimitvType.add(st.nextToken());
				// System.out.println(pPrimitvType);
			}
		}

		// --------------------------
		// List list to Array
		// --------------------------

		//
		// Vertices
		//
		vertices = new float[verticesList.size()][];

		for (int i = 0; i < verticesList.size(); i++) {
			vertices[i] = new float[verticesList.get(i).length];

			vertices[i][0] = verticesList.get(i)[0];
			vertices[i][1] = -verticesList.get(i)[1];
			vertices[i][2] = -verticesList.get(i)[2];
		}
		System.out.println("vertices");
		System.out.println(Arrays.deepToString(vertices));

		verticesMat = new int[verticesMatList.size()];
		for (int i = 0; i < verticesList.size(); i++) {
			verticesMat[i] = verticesMatList.get(i)[0];
		}

		//
		// VertexNormals
		//
		vertexNormals = new float[vertexNormalsList.size()][];
		for (int i = 0; i < vertexNormalsList.size(); i++) {
			vertexNormals[i] = new float[vertexNormalsList.get(i).length];
			for (int j = 0; j < vertexNormalsList.get(i).length; j++) {
				vertexNormals[i][j] = vertexNormalsList.get(i)[j];
			}
		}
		System.out.println("vertexNormals");
		System.out.println(Arrays.deepToString(vertexNormals));

		//
		// faces to triangels
		//
		triangles = new int[facesV.size()][];
		for (int i = 0; i < facesV.size(); i++) {
			triangles[i] = new int[facesV.get(i).length];
			// triangles[i][0] = facesV.get(i)[0] - 1;
			// triangles[i][1] = facesV.get(i)[1] - 1;
			// triangles[i][2] = facesV.get(i)[2] - 1;

			for (int j = 0; j < facesV.get(i).length; j++) {
				triangles[i][j] = facesV.get(i)[j] - 1;
			}
		}
		System.out.println("triangles");
		System.out.println(Arrays.deepToString(triangles));

		//
		// Material
		//
		materials = new float[textureCounter][];
		for (int i = 0; i < textureCounter; i++) {
			materials[i] = new float[] { Ka.get(i)[0], Ka.get(i)[1], Ka.get(i)[2], Kd.get(i)[0], Kd.get(i)[1],
					Kd.get(i)[2], Ks.get(i)[0], Ks.get(i)[1], Ks.get(i)[2] };
		}
		System.out.println("materials");
		System.out.println(Arrays.deepToString(materials));

		materialsN = new int[textureCounter];
		for (int i = 0; i < textureCounter; i++) {
			materialsN[i] = (int) Ns.get(i)[0];
		}
		System.out.println();
		System.out.println(Arrays.toString(materialsN));

		// BBox berechnen
		calcBoundingBox();
	}

	// ===========================================================
	// Method
	// ===========================================================

	/**
	 * Open and read a file, and return the lines in the file as a list of Strings.
	 */
	private List<String> readFile(String filename) {
		List<String> records = new ArrayList<String>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			String line;
			while ((line = reader.readLine()) != null) {
				records.add(line);
			}
			reader.close();
			return records;
		} catch (Exception e) {
			System.err.format("Exception occurred trying to read '%s'.", filename);
			e.printStackTrace();
			return null;
		}
	}
}
