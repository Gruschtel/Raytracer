package de.fhkl.imst.i.cgma.raytracer.importer;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.ShortBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import de.fhkl.imst.i.cgma.raytracer.file.RT_Object;

/**
 * GLTF_Mesh
 */
public class GLTF_Mesh extends RT_Object {

//												  r_a   g_a   b_a   r_a   g_d   b_d   r_a   g_s   b_s
	private float[] basisMaterial = new float[] { 0.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.5f, 0.3f, 0.3f, 0.3f };
	private int basisMaterialN = 10;

	public float[][] materials;
	public int[] materialsN;

	public float[][] vertices;
	public int[] verticesMat;

	public int[][] triangles;

	public char fgp = 'f'; // flat, gouraud, phong

	// calculated information
	public float[][] vertexNormals;
	public float[][] vertexColors;

	public float[][] triangleNormals;
	public float[][] triangleColors;
	public float[] triangleAreas;

	@Override
	public String getHeader() {
		return "GLTF";
	}

	@Override
	public void readContent(File f) throws IOException {
		ObjectMapper mapper = new ObjectMapper();
		JsonNode rootNode = null;

		try {

			rootNode = mapper.readTree(f);
			loadMeshes(rootNode);

		} catch (JsonProcessingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void loadMeshes(JsonNode rootNode) {
		List<ByteBuffer> buffers = new ArrayList<>();
		List<float[]> materialsList = new ArrayList<>();

		List<float[]> tmpVertices = new ArrayList<>();
		List<int[]> tmpVerticesMaterial = new ArrayList<>();
		List<float[]> tmpNormals = new ArrayList<>();
		List<int[]> tmpTriangles = new ArrayList<>();
		int triangleOffset = 0;

//		Get all buffers
		for (JsonNode jsonNode : rootNode.get("buffers")) {
			String uri = jsonNode.get("uri").asText();
			String data = (String) uri.subSequence(uri.lastIndexOf("base64,") + 7, uri.length());
			byte[] bufferData = Base64.getDecoder().decode(data);
			ByteBuffer bb = ByteBuffer.allocateDirect(bufferData.length);
			bb.put(bufferData);
			bb.order(ByteOrder.LITTLE_ENDIAN);
			bb.rewind();
			buffers.add(bb);
		}

//		Get all material
		if (rootNode.has("materials")) {
			for (JsonNode jsonNode : rootNode.get("materials")) {
				JsonNode diffuseNode = jsonNode.get("pbrMetallicRoughness").get("baseColorFactor");
				if (diffuseNode != null) {
					float[] material = new float[9];
					float[] diffuseColor = new float[diffuseNode.size()];
					for (int i = 0; i < diffuseColor.length; i++) {
						diffuseColor[i] = diffuseNode.get(i).floatValue();
					}

					material[0] = basisMaterial[0];
					material[1] = basisMaterial[1];
					material[2] = basisMaterial[2];
					material[3] = diffuseColor[0];
					material[4] = diffuseColor[1];
					material[5] = diffuseColor[2];
					material[6] = basisMaterial[6];
					material[7] = basisMaterial[7];
					material[8] = basisMaterial[8];

					materialsList.add(material);
				} else {
					materialsList.add(basisMaterial);
				}
			}
		} else {
			materialsList.add(basisMaterial);
		}

		int meshIndex = 0;
		for (JsonNode jsonNode : rootNode.get("meshes")) {

//			Get transformation matrix
			float[][] translation = new float[][] { { 1f, 0f, 0f, 0f }, { 0f, 1f, 0f, 0f }, { 0f, 0f, 1f, 0f },
					{ 0f, 0f, 0f, 1f } };
			float[][] rotation = new float[][] { { 1f, 0f, 0f, 0f }, { 0f, 1f, 0f, 0f }, { 0f, 0f, 1f, 0f },
					{ 0f, 0f, 0f, 1f } };
			float[][] scale = new float[][] { { 1f, 0f, 0f, 0f }, { 0f, 1f, 0f, 0f }, { 0f, 0f, 1f, 0f },
					{ 0f, 0f, 0f, 1f } };

//			Rotation
			if (rootNode.get("nodes").get(meshIndex).has("rotation")) {
				JsonNode rotNode = rootNode.get("nodes").get(meshIndex).get("rotation");

				float qx = rotNode.get(0).floatValue();
				float qy = rotNode.get(1).floatValue();
				float qz = rotNode.get(2).floatValue();
				float qw = rotNode.get(3).floatValue();

				float sqx = qx * qx;
				float sqy = qy * qy;
				float sqz = qz * qz;
				float sqw = qw * qw;

				float norm = 1 / (sqx + sqy + sqz + sqw);

				rotation[0][0] = (sqx - sqy - sqz + sqw) * norm;
				rotation[1][1] = (-sqx + sqy - sqz + sqw) * norm;
				rotation[2][2] = (-sqx - sqy + sqz + sqw) * norm;

				float tmp1 = qx * qy;
				float tmp2 = qz * qw;
				rotation[1][0] = 2 * (tmp1 + tmp2) * norm;
				rotation[0][1] = 2 * (tmp1 - tmp2) * norm;

				tmp1 = qx * qz;
				tmp2 = qy * qw;
				rotation[2][0] = 2 * (tmp1 - tmp2) * norm;
				rotation[0][2] = 2 * (tmp1 + tmp2) * norm;

				tmp1 = qy * qz;
				tmp2 = qx * qw;
				rotation[2][1] = 2 * (tmp1 + tmp2) * norm;
				rotation[1][2] = 2 * (tmp1 - tmp2) * norm;

			}

//			Translation
			if (rootNode.get("nodes").get(meshIndex).has("translation")) {
				JsonNode transNode = rootNode.get("nodes").get(meshIndex).get("translation");
				translation[0][3] = transNode.get(0).floatValue();
				translation[1][3] = transNode.get(1).floatValue();
				translation[2][3] = transNode.get(2).floatValue();
			}

//			Scale
			if (rootNode.get("nodes").get(meshIndex).has("scale")) {
				JsonNode scaNode = rootNode.get("nodes").get(meshIndex).get("scale");
				scale[0][0] = scaNode.get(0).floatValue();
				scale[1][1] = scaNode.get(1).floatValue();
				scale[2][2] = scaNode.get(2).floatValue();
			}

			float[][] transformationMatrix = matrixMult(matrixMult(translation, rotation), scale);


			int primitiveOffset = 0;
			
			List<Integer> tmpIndexList = new ArrayList<>();
			
			for (JsonNode primitive : jsonNode.get("primitives")) {

//			Get vertices
				int vertexAccessor = primitive.findValue("attributes").get("POSITION").asInt();

				int vertexBufferView = rootNode.get("accessors").get(vertexAccessor).get("bufferView").asInt();
				int vertexBuffer = rootNode.get("bufferViews").get(vertexBufferView).get("buffer").asInt();
				int vertexByteOffset = rootNode.get("bufferViews").get(vertexBufferView).get("byteOffset").asInt();
				int vertexByteLength = rootNode.get("bufferViews").get(vertexBufferView).get("byteLength").asInt();

				FloatBuffer fb = buffers.get(vertexBuffer).asFloatBuffer();
				fb.position(vertexByteOffset / Float.BYTES);
				fb.limit((vertexByteOffset + vertexByteLength) / Float.BYTES);

				float[] floatArray = new float[fb.remaining()];
				fb.get(floatArray, 0, floatArray.length);

//			Transform all vertices
				for (int i = 0; i < floatArray.length; i += 3) {
					float[] vertex = transformVector(transformationMatrix,
							new float[] { floatArray[i], floatArray[i + 1], floatArray[i + 2] });
					floatArray[i] = vertex[0];
					floatArray[i + 1] = vertex[1];
					floatArray[i + 2] = vertex[2];

				}
				
				tmpVertices.add(floatArray);

//			Get triangles
				int triangleAccessor = primitive.findValue("indices").asInt();

				int triangleBufferView = rootNode.get("accessors").get(triangleAccessor).get("bufferView").asInt();
				int triangleBuffer = rootNode.get("bufferViews").get(triangleBufferView).get("buffer").asInt();
				int triangleByteOffset = rootNode.get("bufferViews").get(triangleBufferView).get("byteOffset").asInt();
				int triangleByteLength = rootNode.get("bufferViews").get(triangleBufferView).get("byteLength").asInt();

				ShortBuffer sb = buffers.get(triangleBuffer).asShortBuffer();
				sb.position(triangleByteOffset / Short.BYTES);
				sb.limit((triangleByteOffset + triangleByteLength) / Short.BYTES);

				short[] shortArray = new short[sb.remaining()];
				sb.get(shortArray, 0, shortArray.length);
				
				int[] outArray = new int[shortArray.length];
				for (int i = 0; i < shortArray.length; i++) {
					outArray[i] = shortArray[i];

					outArray[i] += triangleOffset + primitiveOffset;
					tmpIndexList.add(outArray[i]);
				}

//				Get Material
				int[] vertMats = new int[floatArray.length / 3];
				JsonNode material = primitive.findValue("material");
				if (material != null) { // Check if there is a material
					int materialIndex = material.asInt();
					for (int i = 0; i < vertMats.length; i++) {
						vertMats[i] = materialIndex;
					}
				} else {
					for (int i = 0; i < vertMats.length; i++) {
						vertMats[i] = 0;
					}
				}
				
				tmpVerticesMaterial.add(vertMats);

				int max = 0;
				for (int is : outArray) {
					max = Math.max(is, max);
				}

				primitiveOffset = max - triangleOffset + 1 ;
			}
			
			int max = 0;
			int[] tOut = new int[tmpIndexList.size()];
			for (int i = 0; i < tmpIndexList.size(); i++) {
				tOut[i] = tmpIndexList.get(i);
			}
			tmpTriangles.add(tOut);
			for (int is : tOut) {
				max = Math.max(is, max);
			}
			
			triangleOffset = max + 1;

			meshIndex++;
		}

//		Combine the arrays
		materials = new float[materialsList.size()][];
		materialsN = new int[materialsList.size()];
		for (int i = 0; i < materialsList.size(); i++) {
			materials[i] = materialsList.get(i);
			materialsN[i] = basisMaterialN;
		}

		ArrayList<Integer> vmList = new ArrayList<>();
		for (int[] tvm : tmpVerticesMaterial) {
			for (int i : tvm) {
				vmList.add(i);
			}
		}
		verticesMat = ArrayUtils.toPrimitive(vmList.toArray(new Integer[vmList.size()]));

		ArrayList<Float> vList = new ArrayList<>();
		for (float[] tv : tmpVertices) {
			for (float f : tv) {
				vList.add(f);
			}
		}
		vertices = new float[vList.size() / 3][];
		for (int i = 0, j = 0; i < vertices.length; i++, j += 3) {
			vertices[i] = new float[] { vList.get(j), vList.get(j + 1), vList.get(j + 2) };
		}

		ArrayList<Integer> tList = new ArrayList<>();
		for (int[] tt : tmpTriangles) {
			for (int i : tt) {
				tList.add(i);
			}
		}
		triangles = new int[tList.size() / 3][];
		for (int i = 0, j = 0; i < triangles.length; i++, j += 3) {
			triangles[i] = new int[] { tList.get(j), tList.get(j + 1), tList.get(j + 2) };
		}

		calcBoundingBox();
	}

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
	}

	private float[][] matrixMult(float[][] matA, float[][] matB) {
		int aRows = matA.length;
		int aColumns = matA[0].length;
		int bRows = matB.length;
		int bColumns = matB[0].length;

		if (aColumns != bRows) {
			throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
		}

		float[][] C = new float[aRows][bColumns];
		for (int i = 0; i < aRows; i++) {
			for (int j = 0; j < bColumns; j++) {
				for (int k = 0; k < aColumns; k++) { // aColumn
					C[i][j] += matA[i][k] * matB[k][j];
				}
			}
		}

		return C;
	}

	private float[] transformVector(float[][] matrix, float[] vector) {
		float[] inVector = new float[] { vector[0], vector[1], vector[2], 1 };

		int matRows = matrix.length;
		int matCols = matrix[0].length;
		if (inVector.length != matCols)
			throw new RuntimeException("Illegal matrix dimensions.");

		float[] y = new float[matRows];
		for (int i = 0; i < matRows; i++)
			for (int j = 0; j < matCols; j++)
				y[i] += matrix[i][j] * inVector[j];
		return y;
	}

	@Override
	protected void readContent(LineNumberReader f) throws IOException {

	}
}
