package de.fhkl.imst.i.cgma.raytracer;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;

import de.fhkl.imst.i.cgma.raytracer.file.GLTF_Mesh;
import de.fhkl.imst.i.cgma.raytracer.file.I_Sphere;
import de.fhkl.imst.i.cgma.raytracer.file.Obj_Mesh;
import de.fhkl.imst.i.cgma.raytracer.file.RTFile;
import de.fhkl.imst.i.cgma.raytracer.file.RTFileReader;
import de.fhkl.imst.i.cgma.raytracer.file.RT_Object;
import de.fhkl.imst.i.cgma.raytracer.file.T_Mesh;
import de.fhkl.imst.i.cgma.raytracer.gui.IRayTracerImplementation;
import de.fhkl.imst.i.cgma.raytracer.gui.RayTracerGui;

public class Raytracer00 implements IRayTracerImplementation {
	// viewing volume with infinite end
	private float fovyDegree;
	private float near;
	private float fovyRadians;

	// one hardcoded point light as a minimal solution :-(
	private float[] Ia = { 0.25f, 0.25f, 0.25f }; // ambient light color
	private float[] Ids = { 1.0f, 1.0f, 1.0f }; // diffuse and specular light
	// color
	private float[] ICenter = { 4.0f, 4.0f, 2.0f }; // center of point light

	RayTracerGui gui = new RayTracerGui(this);

	private int resx, resy; // viewport resolution
	private float h, w, aspect; // window height, width and aspect ratio

	Vector<RT_Object> objects;

	private Raytracer00() {
		try {
			String directory = System.getProperty("user.dir");
//		    gui.addObject(RTFileReader.read(I_Sphere.class, new File("data/ikugel.dat")));
//			gui.addObject(RTFileReader.read(T_Mesh.class, new File("data/dreieck1.dat")));

//			gui.addObject(RTFileReader.read(Obj_Mesh.class, new File(directory + "/data/object2.obj")));
			gui.addObject(RTFileReader.read(Obj_Mesh.class, new File(directory + "/data/obj_moreMeshes.obj")));
			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/obj_elefant.gltf")));

//			gui.addObject(RTFileReader.read(Obj_Mesh.class, new File(directory + "/data/sphere1.obj")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/simpleTriangle.gltf")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/untitled.gltf")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/sphere.gltf")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/sphere1.gltf")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/elefant1.gltf")));
//			gui.addObject(RTFileReader.read(GLTF_Mesh.class, new File(directory + "/data/Bender.gltf")));

			objects = gui.getObjects();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void setViewParameters(float fovyDegree, float near) {
		// set attributes fovyDegree, fovyRadians, near
		this.fovyDegree = fovyDegree;
		this.fovyRadians = (float) (this.fovyDegree * Math.PI / 180f);
		this.near = near;

		// set attributes resx, resy, aspect
		resx = gui.getResX();
		resy = gui.getResY();
		aspect = resx / (float) resy;

		// set attributes h, w
		h = (float) (2 * near * Math.tan(fovyRadians / 2f));
		w = h * aspect;

//		System.out.println("Aspect: " + aspect);
//		System.out.println("fovyradians: " + fovyRadians);
//		System.out.println("h: " + h);
//		System.out.println("w: " + w);
	}

	@Override
	public void doRayTrace() {
		float x, y, z; // intersection point in viewing plane
		float rayEx, rayEy, rayEz; // eye point==ray starting point
		float rayVx, rayVy, rayVz; // ray vector

		// prepare mesh data (normals and areas)
		prepareMeshData();

		// hardcoded viewing volume with fovy and near
		setViewParameters(90.0f, 1.0f);
		// set eye point
		rayEx = 0.0f;
		rayEy = 0.0f;
		rayEz = 0.0f;

		z = -near;

		Color color;

		// prepare mesh data for shading
		precalculateMeshDataShading();

		Random rd = new Random();
		// xp, yp: pixel coordinates
		for (int xp = 0; xp < resx; ++xp) {
			for (int yp = 0; yp < resy; ++yp) {
				gui.updateSTime();
				// for demo purposes
				gui.setPixel(xp, yp, Color.WHITE.getRGB());
				// gui.setPixel(xp, yp, new Color(rd.nextFloat(), rd.nextFloat(),
				// rd.nextFloat()).getRGB());

				// x, y: view coordinates

				x = xp / (resx - 1f) * w - w / 2f;
				y = (resy - yp - 1f) / (resy - 1f) * h - h / 2f;
				// ray vector
				rayVx = x - rayEx;
				rayVy = y - rayEy;
				rayVz = z - rayEz;
				// get color or null along the ray
				// color=traceRayAndGetColor...
				// if color!=null set pixel with color
				color = traceRayAndGetColor(rayEx, rayEy, rayEz, rayVx, rayVy, rayVz);
				if (color != null) {
					gui.setPixel(xp, yp, color.getRGB());
					// System.out.println("gui.setPixel" + color.getRed());
				}
			}
		}
		System.out.println("Render Finished!");
	}

	// returns Color object or null if no intersection was found
	private Color traceRayAndGetColor(float rayEx, float rayEy, float rayEz, float rayVx, float rayVy, float rayVz) {
		// RTFile scene = gui.getFile();

		double minT = Float.MAX_VALUE;
		int minObjectsIndex = -1;
		int minIndex = -1;
		float[] minIP = new float[3];
		float[] minN = new float[3];
		float[] minMaterial = new float[3];
		float minMaterialN = 1;
		float bu = 0, bv = 0, bw = 1;

		float[] v = new float[3];
		float[] l = new float[3];

		RTFile scene;
		I_Sphere sphere;
		T_Mesh mesh;
		Obj_Mesh obj_mesh;
		GLTF_Mesh gltf_Mesh;

		float[] rv = { rayVx, rayVy, rayVz };
		float[] re = { rayEx, rayEy, rayEz };
		float[] m;
		float[] me = new float[3];

		// viewing vector at intersection point
		v[0] = -rv[0];
		v[1] = -rv[1];
		v[2] = -rv[2];

		// loop over all scene objects to find the nearest intersection, that
		// is:
		// object with number minObjectIndex
		// minT is the minimal factor t of the ray equation s(t)=rayE+t*rayV
		// where the nearest intersection takes place
		for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
			scene = objects.get(objectsNumber);

			// object is an implicit sphere?
			if (scene instanceof I_Sphere) {
				sphere = (I_Sphere) scene;

				m = sphere.center;
				me[0] = re[0] - m[0];
				me[1] = re[1] - m[1];
				me[2] = re[2] - m[2];

				float t;

				// no bounding box hit? -> next object
				if (!bboxHit(sphere, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz))
					continue;

				// ray intersection uses quadratic equation
				float a, b, c, d;

				a = rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2];

				b = 2 * (rv[0] * (me[0]) + rv[1] * (me[1]) + rv[2] * (me[2]));

				c = (me[0] * me[0] + me[1] * me[1] + me[2] * me[2]) - sphere.radius * sphere.radius;

				// positive discriminant determines intersection b^2 - 4ac
				d = b * b - 4 * a * c;

				// no intersection point? => next object
				if (d <= 0)
					continue;
				// from here: intersection takes place!

				// calculate first intersection point with sphere along the
				// ray

				t = (float) ((-b - Math.sqrt(d)) / (2 * a));
				// already a closer intersection point? => next object
				if (t >= minT)
					continue;

				// from here: t < minT

				// I'm the winner until now!

				minT = t;
				minObjectsIndex = objectsNumber;

				// prepare everything for phong shading

				// the intersection point
				minIP[0] = rayEx + t * rayVx;
				minIP[1] = rayEy + t * rayVy;
				minIP[2] = rayEz + t * rayVz;

				// the normal vector at the intersection point
				minN[0] = minIP[0] - sphere.center[0];
				minN[1] = minIP[1] - sphere.center[1];
				minN[2] = minIP[2] - sphere.center[2];

				// the material
				minMaterial = sphere.material;
				minMaterialN = sphere.materialN;

			}

			// object is a triangle mesh?
			else if (scene instanceof T_Mesh) {
				mesh = (T_Mesh) scene;

				float t;
				float[] n;
				float[] ip = new float[3];

				// no bounding box hit? -> next object
				if (!bboxHit(mesh, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz))
					continue;

				float a, rayVn, pen;
				float[] p1, p2, p3;
				float[] ai = new float[3];

				// loop over all triangles
				for (int i = 0; i < mesh.triangles.length; i++) {
					// get the three vertices
					p1 = mesh.vertices[mesh.triangles[i][0]];
					p2 = mesh.vertices[mesh.triangles[i][1]];
					p3 = mesh.vertices[mesh.triangles[i][2]];
//
//					// intermediate version
//					// calculate normal n and triangle area a
//					n = new float[3];
//					a = calculateN(n, p1, p2, p3);

					// fetch precalculated face areas and face normals
					a = mesh.triangleAreas[i];
					n = mesh.triangleNormals[i];

					rayVn = rv[0] * n[0] + rv[1] * n[1] + rv[2] * n[2];

					// backface? => next triangle
					if (rayVn >= 0)
						continue;
					// no intersection point? => next triangle
					if (Math.abs(rayVn) < 1E-7)
						continue;

					pen = (p1[0] - rayEx) * n[0] + (p1[1] - rayEy) * n[1] + (p1[2] - rayEz) * n[2];

					// calculate intersection point with plane along the ray
					t = pen / rayVn;

					// already a closer intersection point? => next triangle
					if (t >= minT)
						continue;

					// the intersection point with the plane
					ip[0] = rayEx + (float) t * rayVx;
					ip[1] = rayEy + (float) t * rayVy;
					ip[2] = rayEz + (float) t * rayVz;

					// no intersection point with the triangle? => next
					// triangle
					if (!triangleTest(ip, p1, p2, p3, a, ai))
						continue;

					// from here: t < minT and triangle intersection
					// I'm the winner until now!

					minT = t;
					minObjectsIndex = objectsNumber;
					minIndex = i;

					// prepare everything for shading alternatives

					// the intersection point
					minIP[0] = ip[0];
					minIP[1] = ip[1];
					minIP[2] = ip[2];

					switch (mesh.fgp) {
					case 'f':
					case 'F':

//						// the normal is the surface normal
//						minN[0] = n[0];
//						minN[1] = n[1];
//						minN[2] = n[2];
//
//						// the material is the material of the first triangle point
//						int matIndex = mesh.verticesMat[mesh.triangles[minIndex][0]];
//						minMaterial = mesh.materials[matIndex];
//						minMaterialN = mesh.materialsN[matIndex];

						break;
					case 'g':
					case 'G':
						// remember barycentric coordinates bu, bv, bw for shading
						bu = ai[0];
						bv = ai[1];
						bw = ai[2];

						break;
					case 'p':
					case 'P':
						// the normal is barycentrically interpolated between
						// the three vertices
						bu = ai[0];
						bv = ai[1];
						bw = ai[2];

						float nTemp[] = new float[3];
						nTemp[0] = mesh.vertexNormals[mesh.triangles[minIndex][2]][0] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][0] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][0] * bw;
						nTemp[1] = mesh.vertexNormals[mesh.triangles[minIndex][2]][1] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][1] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][1] * bw;
						nTemp[2] = mesh.vertexNormals[mesh.triangles[minIndex][2]][2] * bu
								+ mesh.vertexNormals[mesh.triangles[minIndex][0]][2] * bv
								+ mesh.vertexNormals[mesh.triangles[minIndex][1]][2] * bw;

						minN = nTemp;

						// intermediate version
						// the material is not interpolated
						// matIndex =
						// mesh.verticesMat[mesh.triangles[minIndex][0]];
						// minMaterial = mesh.materials[matIndex];
						// minMaterialN = mesh.materialsN[matIndex];

						// the material is barycentrically interpolated between
						// the three vertex materials
						int matIndex0 = mesh.verticesMat[mesh.triangles[minIndex][0]];
						int matIndex1 = mesh.verticesMat[mesh.triangles[minIndex][1]];
						int matIndex2 = mesh.verticesMat[mesh.triangles[minIndex][2]];
						float materialTemp[] = new float[9];
						int materialNTemp;
						for (int k = 0; k < 9; k++) {
							materialTemp[k] = mesh.materials[matIndex0][k] * bu + mesh.materials[matIndex1][k] * bv
									+ mesh.materials[matIndex2][k] * bw;
						}

						minMaterial = materialTemp;
						materialNTemp = (int) (mesh.materialsN[matIndex0] * bu + mesh.materialsN[matIndex1] * bv
								+ mesh.materialsN[matIndex2] * bw);
						minMaterialN = materialNTemp;
					}
				}
			} else if (scene instanceof Obj_Mesh) {
				// TODO: Obj_Mesh 1
				obj_mesh = (Obj_Mesh) scene;

				float t;
				float[] n;
				float[] ip = new float[3];

				// no bounding box hit? -> next object
				if (!bboxHit(obj_mesh, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz))
					continue;


				float a, rayVn, pen;
				float[] p1, p2, p3;
				float[] ai = new float[3];

				for (int i = 0; i < obj_mesh.triangles.length; i++) {
					// get the three vertices
					p1 = obj_mesh.vertices[obj_mesh.triangles[i][0]];
					p2 = obj_mesh.vertices[obj_mesh.triangles[i][1]];
					p3 = obj_mesh.vertices[obj_mesh.triangles[i][2]];
//
//					// intermediate version
//					// calculate normal n and triangle area a
//					n = new float[3];
//					a = calculateN(n, p1, p2, p3);

					// fetch precalculated face areas and face normals
					a = obj_mesh.triangleAreas[i];
					n = obj_mesh.triangleNormals[i];

					rayVn = rv[0] * n[0] + rv[1] * n[1] + rv[2] * n[2];

					// backface? => next triangle
					if (rayVn >= 0)
						continue;
					// no intersection point? => next triangle
					if (Math.abs(rayVn) < 1E-7)
						continue;

					pen = (p1[0] - rayEx) * n[0] + (p1[1] - rayEy) * n[1] + (p1[2] - rayEz) * n[2];

					// calculate intersection point with plane along the ray
					t = pen / rayVn;

					// already a closer intersection point? => next triangle
					if (t >= minT)
						continue;

					// the intersection point with the plane
					ip[0] = rayEx + (float) t * rayVx;
					ip[1] = rayEy + (float) t * rayVy;
					ip[2] = rayEz + (float) t * rayVz;

					// no intersection point with the triangle? => next
					// triangle
					if (!triangleTest(ip, p1, p2, p3, a, ai))
						continue;

					// from here: t < minT and triangle intersection
					// I'm the winner until now!

					minT = t;
					minObjectsIndex = objectsNumber;
					minIndex = i;

					// prepare everything for shading alternatives

					// the intersection point
					minIP[0] = ip[0];
					minIP[1] = ip[1];
					minIP[2] = ip[2];

					switch (obj_mesh.fgp) {
					case 'f':
					case 'F':

//						// the normal is the surface normal
//						minN[0] = n[0];
//						minN[1] = n[1];
//						minN[2] = n[2];
//
//						// the material is the material of the first triangle point
//						int matIndex = obj_mesh.verticesMat[obj_mesh.triangles[minIndex][0]];
//						minMaterial = obj_mesh.materials[matIndex];
//						minMaterialN = obj_mesh.materialsN[matIndex];

						break;
					}
				}

			} else if (scene instanceof GLTF_Mesh) {
				gltf_Mesh = (GLTF_Mesh) scene;

				float t;
				float[] n;
				float[] ip = new float[3];

				// no bounding box hit? -> next object
				if (!bboxHit(gltf_Mesh, rayEx, rayEy, rayEz, rayVx, rayVy, rayVz))
					continue;


				float a, rayVn, pen;
				float[] p1, p2, p3;
				float[] ai = new float[3];

				for (int i = 0; i < gltf_Mesh.triangles.length; i++) {
					// get the three vertices
					p1 = gltf_Mesh.vertices[gltf_Mesh.triangles[i][0]];
					p2 = gltf_Mesh.vertices[gltf_Mesh.triangles[i][1]];
					p3 = gltf_Mesh.vertices[gltf_Mesh.triangles[i][2]];
//
//					// intermediate version
//					// calculate normal n and triangle area a
//					n = new float[3];
//					a = calculateN(n, p1, p2, p3);

					// fetch precalculated face areas and face normals
					a = gltf_Mesh.triangleAreas[i];
					n = gltf_Mesh.triangleNormals[i];

					rayVn = rv[0] * n[0] + rv[1] * n[1] + rv[2] * n[2];

					// backface? => next triangle
					if (rayVn >= 0)
						continue;
					// no intersection point? => next triangle
					if (Math.abs(rayVn) < 1E-7)
						continue;

					pen = (p1[0] - rayEx) * n[0] + (p1[1] - rayEy) * n[1] + (p1[2] - rayEz) * n[2];

					// calculate intersection point with plane along the ray
					t = pen / rayVn;

					// already a closer intersection point? => next triangle
					if (t >= minT)
						continue;

					// the intersection point with the plane
					ip[0] = rayEx + (float) t * rayVx;
					ip[1] = rayEy + (float) t * rayVy;
					ip[2] = rayEz + (float) t * rayVz;

					// no intersection point with the triangle? => next
					// triangle
					if (!triangleTest(ip, p1, p2, p3, a, ai))
						continue;

					// from here: t < minT and triangle intersection
					// I'm the winner until now!

					minT = t;
					minObjectsIndex = objectsNumber;
					minIndex = i;

					// prepare everything for shading alternatives

					// the intersection point
					minIP[0] = ip[0];
					minIP[1] = ip[1];
					minIP[2] = ip[2];

					switch (gltf_Mesh.fgp) {
					case 'f':
					case 'F':

//						// the normal is the surface normal
//						minN[0] = n[0];
//						minN[1] = n[1];
//						minN[2] = n[2];
//
//						// the material is the material of the first triangle point
//						int matIndex = obj_mesh.verticesMat[obj_mesh.triangles[minIndex][0]];
//						minMaterial = obj_mesh.materials[matIndex];
//						minMaterialN = obj_mesh.materialsN[matIndex];

						break;
					}
				}
			} else
				continue; // return null;
		}

		// no intersection point found => return with no result
		if (minObjectsIndex == -1)
			return null;

		// light vector at the intersection point
		l[0] = ICenter[0] - minIP[0];
		l[1] = ICenter[1] - minIP[1];
		l[2] = ICenter[2] - minIP[2];

		normalize(minN);
		normalize(v);
		normalize(l);

		// decide which shading model will be applied

		// implicit: only phong shading available => shade=illuminate
		if (objects.get(minObjectsIndex) instanceof I_Sphere) {
			return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);
		}
		// triangle mesh: flat, gouraud or phong shading according to file data
		else if (objects.get(minObjectsIndex).getHeader() == "TRIANGLE_MESH") {
			mesh = ((T_Mesh) objects.get(minObjectsIndex));
			switch (mesh.fgp) {
			case 'f':
			case 'F':
//			// illumination can be calculated here
//			// this is a variant between flat und phong shading
//			return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);

				// lookup triangle color of triangle hit
				return new Color(mesh.triangleColors[minIndex][0], mesh.triangleColors[minIndex][1],
						mesh.triangleColors[minIndex][2]);
			case 'g':
			case 'G':
				// the color is barycentrically interpolated between the three
				// vertex colors
				float colorf[] = new float[3];
				colorf[0] = mesh.vertexColors[mesh.triangles[minIndex][2]][0] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][0] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][0] * bw;
				colorf[1] = mesh.vertexColors[mesh.triangles[minIndex][2]][1] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][1] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][1] * bw;
				colorf[2] = mesh.vertexColors[mesh.triangles[minIndex][2]][2] * bu
						+ mesh.vertexColors[mesh.triangles[minIndex][0]][2] * bv
						+ mesh.vertexColors[mesh.triangles[minIndex][1]][2] * bw;

				return new Color(colorf[0] > 1 ? 1 : colorf[0], colorf[1] > 1 ? 1 : colorf[1],
						colorf[2] > 1 ? 1 : colorf[2]);
			case 'p':
			case 'P':
				// calculate the color per per pixel phong lightning
				return phongIlluminate(minMaterial, minMaterialN, l, minN, v, Ia, Ids);
			// return new Color(material[3], material[4], material[5]);
			// break;

			}

		} else if (objects.get(minObjectsIndex).getHeader() == "Obj") {
			// TODO: Obj_Mesh 2
			obj_mesh = ((Obj_Mesh) objects.get(minObjectsIndex));

			switch (obj_mesh.fgp) {
			case 'f':
			case 'F':

				return new Color(obj_mesh.triangleColors[minIndex][0], obj_mesh.triangleColors[minIndex][1],
						obj_mesh.triangleColors[minIndex][2]);
			}
		} else if (objects.get(minObjectsIndex).getHeader() == "GLTF") {

			gltf_Mesh = ((GLTF_Mesh) objects.get(minObjectsIndex));

			switch (gltf_Mesh.fgp) {
			case 'f':
			case 'F':

				return new Color(gltf_Mesh.triangleColors[minIndex][0], gltf_Mesh.triangleColors[minIndex][1],
						gltf_Mesh.triangleColors[minIndex][2]);

			default:
				break;
			}
		}
		return null;
	}

	// calculate phong illumination model with material parameters material and
	// materialN, light vector l, normal vector n, viewing vector v, ambient
	// light Ia, diffuse and specular light Ids
	// return value is a new Color object
	private Color phongIlluminate(float[] material, float materialN, float[] l, float[] n, float[] v, float[] Ia,
			float[] Ids) {
		float ir = 0, ig = 0, ib = 0; // reflected intensity, rgb channels
		float[] r = new float[3]; // reflection vector
		float ln, rv; // scalar products <l,n> and <r,v>

		// <l,n>

		ln = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];

		// ambient component, Ia*ra
		ir += Ia[0] * material[0];
		ig += Ia[1] * material[1];
		ib += Ia[2] * material[2];

		// diffuse component, Ids*rd*<l,n>
		if (ln > 0) {
			ir += Ids[0] * material[3] * ln;
			ig += Ids[1] * material[4] * ln;
			ib += Ids[2] * material[5] * ln;

			// reflection vector r=2*<l,n>*n-l
			r[0] = 2 * ln * n[0] - l[0];
			r[1] = 2 * ln * n[1] - l[1];
			r[2] = 2 * ln * n[2] - l[2];

			normalize(r);
			// <r,v>
			rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

			// specular component, Ids*rs*<r,v>^n
			if (rv > 0) {
				float pow = (float) Math.pow(rv, materialN);
				ir += Ids[0] * material[6] * pow;
				ig += Ids[1] * material[7] * pow;
				ib += Ids[2] * material[8] * pow;
			}
		}

		return new Color(ir > 1 ? 1 : ir, ig > 1 ? 1 : ig, ib > 1 ? 1 : ib);
	}

	public float normalize(float[] vector3D) {
		float length = (float) Math
				.sqrt(vector3D[0] * vector3D[0] + vector3D[1] * vector3D[1] + vector3D[2] * vector3D[2]);

		vector3D[0] = vector3D[0] / length;
		vector3D[1] = vector3D[1] / length;
		vector3D[2] = vector3D[2] / length;

		return length;
	}

	// calculate normalized face normal fn of the triangle p1, p2 and p3
	// the return value is the area of triangle
	// CAUTION: fn is an output parameter; the referenced object will be
	// altered!
	private float calculateN(float[] fn, float[] p1, float[] p2, float[] p3) {
		float ax, ay, az, bx, by, bz;

		// a = Vi2-Vi1, b = Vi3-Vi1
		ax = p2[0] - p1[0];
		ay = p2[1] - p1[1];
		az = p2[2] - p1[2];

		bx = p3[0] - p1[0];
		by = p3[1] - p1[1];
		bz = p3[2] - p1[2];

		// n =a x b
		fn[0] = ay * bz - az * by;
		fn[1] = az * bx - ax * bz;
		fn[2] = ax * by - ay * bx;

		// normalize n, calculate and return area of triangle
		return normalize(fn) / 2f;
	}

	// calculate triangle test
	// is p (the intersection point with the plane through p1, p2 and p3) inside
	// the triangle p1, p2 and p3?
	// the return value answers this question
	// a is an input parameter - the given area of the triangle p1, p2 and p3
	// ai will be computed to be the areas of the sub-triangles to allow to
	// compute barycentric coordinates of the intersection point p
	// ai[0] is associated with bu (p1p2p) across from p3
	// ai[1] is associated with bv (pp2p3) across from p1
	// ai[2] is associated with bw (p1pp3) across form p2
	// CAUTION: ai is an output parameter; the referenced object will be
	// altered!
	private boolean triangleTest(float[] p, float[] p1, float[] p2, float[] p3, float a, float ai[]) {
		float tmp[] = new float[3];

		ai[0] = calculateN(tmp, p1, p2, p) / a;
		ai[1] = calculateN(tmp, p, p2, p3) / a;
		ai[2] = calculateN(tmp, p1, p, p3) / a;


		if (ai[0] + ai[1] + ai[2] <= (1 + 1e-5) && ai[0] + ai[1] + ai[2] >= (1 - 1e-5))
			return true;
		return false;
	}

	// calculate bounding box test
	// decides whether the ray s(t)=rayE+t*rayV intersects the axis aligned
	// t=(s(t)-rayE)/rayV
	// bounding box of object -> return value true
	// six plane intersections with rectangle inside tests; if one succeeds
	// bounding box is hit
	private boolean bboxHit(RT_Object object, float rayEx, float rayEy, float rayEz, float rayVx, float rayVy,
			float rayVz) {
		float t;
		float ip[] = new float[3];

		// front and back
		if (Math.abs(rayVz) > 1E-5) {
			// front xy
			t = (object.min[2] - rayEz) / rayVz;

			ip[0] = rayEx + t * rayVx;
			ip[1] = rayEy + t * rayVy;

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1])
				return true;

			// back xy
			t = (object.max[2] - rayEz) / rayVz;

			ip[0] = rayEx + t * rayVx;
			ip[1] = rayEy + t * rayVy;

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[1] > object.min[1] && ip[1] < object.max[1])
				return true;
		}

		// left and right
		if (Math.abs(rayVx) > 1E-5) {
			// left yz
			t = (object.min[0] - rayEx) / rayVx;

			ip[1] = rayEy + t * rayVy;
			ip[2] = rayEz + t * rayVz;

			if (ip[1] > object.min[1] && ip[1] < object.max[1] && ip[2] > object.min[2] && ip[2] < object.max[2])
				return true;

			// right yz
			t = (object.max[0] - rayEx) / rayVx;

			ip[1] = rayEy + t * rayVy;
			ip[2] = rayEz + t * rayVz;

			if (ip[1] > object.min[1] && ip[1] < object.max[1] && ip[2] > object.min[2] && ip[2] < object.max[2])
				return true;
		}
		// top and bottom
		if (Math.abs(rayVy) > 1E-5) {
			// top xz
			t = (object.min[1] - rayEy) / rayVy;

			ip[0] = rayEx + t * rayVx;
			ip[2] = rayEz + t * rayVz;

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[2] > object.min[2] && ip[2] < object.max[2])
				return true;

			// bottom xz
			t = (object.max[1] - rayEy) / rayVy;

			ip[0] = rayEx + t * rayVx;
			ip[2] = rayEz + t * rayVz;

			if (ip[0] > object.min[0] && ip[0] < object.max[0] && ip[2] > object.min[2] && ip[2] < object.max[2])
				return true;
		}

		return false;
	}

	// precalulation of triangle normals and triangle areas
	private void prepareMeshData() {
		RTFile scene;

		System.out.println("Vorverarbeitung 1 läuft");

		float[] p1, p2, p3;

		for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
			scene = objects.get(objectsNumber);

			if (scene.getHeader() == "TRIANGLE_MESH") {
				T_Mesh mesh = (T_Mesh) scene;

				// init memory
				mesh.triangleNormals = new float[mesh.triangles.length][3];
				mesh.triangleAreas = new float[mesh.triangles.length];

				for (int i = 0; i < mesh.triangles.length; i++) {
					p1 = mesh.vertices[mesh.triangles[i][0]];
					p2 = mesh.vertices[mesh.triangles[i][1]];
					p3 = mesh.vertices[mesh.triangles[i][2]];

					// calculate and store triangle normal n and triangle area a
					mesh.triangleAreas[i] = calculateN(mesh.triangleNormals[i], p1, p2, p3);
				}
			} else if (scene.getHeader() == "Obj") {
				Obj_Mesh mesh = (Obj_Mesh) scene;

				// init memory
				mesh.triangleNormals = new float[mesh.triangles.length][3];
				mesh.triangleAreas = new float[mesh.triangles.length];

				for (int i = 0; i < mesh.triangles.length; i++) {
					p1 = mesh.vertices[(mesh.triangles[i][0])];
					p2 = mesh.vertices[(mesh.triangles[i][1])];
					p3 = mesh.vertices[(mesh.triangles[i][2])];

					// calculate and store triangle normal n and triangle area a
					mesh.triangleAreas[i] = calculateN(mesh.triangleNormals[i], p1, p2, p3);
				}
			} else if (scene.getHeader() == "GLTF") {

				GLTF_Mesh mesh = (GLTF_Mesh) scene;
				mesh.triangleNormals = new float[mesh.triangles.length][3];
				mesh.triangleAreas = new float[mesh.triangles.length];

				for (int i = 0; i < mesh.triangles.length; i++) {
					p1 = mesh.vertices[(mesh.triangles[i][0])];
					p2 = mesh.vertices[(mesh.triangles[i][1])];
					p3 = mesh.vertices[(mesh.triangles[i][2])];

					// calculate and store triangle normal n and triangle area a
					mesh.triangleAreas[i] = calculateN(mesh.triangleNormals[i], p1, p2, p3);
				}
			}
		}
		System.out.println("Vorverarbeitung 1 beendet");
	}

	// view dependend precalculation dependend on type of mesh shading
	// vertexNormals for phong and gouraud shading
	// vertexColors for gouraud shading
	// triangleColors for flat lighting
	private void precalculateMeshDataShading() {

		RTFile scene;

		System.out.println("Vorverarbeitung 2 läuft");

		float rayEx, rayEy, rayEz, rayVx, rayVy, rayVz;
		double rayVn;
		Color color;
		float x, y, z;
		float[] ip = new float[3];
		float[] n = new float[3];
		float[] l = new float[3];
		float[] v = new float[3];
		float[] material;
		float materialN;
		int matIndex;

		for (int objectsNumber = 0; objectsNumber < objects.size(); objectsNumber++) {
			scene = objects.get(objectsNumber);

			if (scene.getHeader() == "TRIANGLE_MESH") {
				T_Mesh mesh = (T_Mesh) scene;

				switch (mesh.fgp) {
				case 'f':
				case 'F':
					// for flat-shading: initialize and calculate triangle
					// colors
					mesh.triangleColors = new float[mesh.triangles.length][3];

					rayEx = 0.0f;
					rayEy = 0.0f;
					rayEz = 0.0f;

					// loop over all triangles
					for (int i = 0; i < mesh.triangles.length; i++) {
						// the intersection point is the first vertex of the
						// triangle
						ip = mesh.vertices[mesh.triangles[i][0]];

						// the material is the material of the first triangle
						// point
						matIndex = mesh.verticesMat[mesh.triangles[i][0]];
						material = mesh.materials[matIndex];
						materialN = mesh.materialsN[matIndex];

						// x, y, z: view coordinates are intersection point
						x = ip[0];
						y = ip[1];
						z = ip[2];

						// ray vector
						rayVx = x - rayEx;
						rayVy = y - rayEy;
						rayVz = z - rayEz;

						// fetch precalculated face normal
						n = mesh.triangleNormals[i];

						rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

						// backface? => next triangle
						if (rayVn >= 0)
							continue;

						// light vector at the intersection point
						l[0] = ICenter[0] - ip[0];
						l[1] = ICenter[1] - ip[1];
						l[2] = ICenter[2] - ip[2];
						normalize(l);

						// viewing vector at intersection point
						v[0] = -rayVx;
						v[1] = -rayVy;
						v[2] = -rayVz;
						normalize(v);

						// illuminate
						color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

						// write color to triangle
						mesh.triangleColors[i][0] = color.getRed() / 255f;
						mesh.triangleColors[i][1] = color.getGreen() / 255f;
						mesh.triangleColors[i][2] = color.getBlue() / 255f;

					}
					break;
				case 'p':
				case 'P':
				case 'g':
				case 'G':
					// initialize and calculate averaged vertex normals
					mesh.vertexNormals = new float[mesh.vertices.length][3];

					// loop over all vertices to initialize
					for (int j = 0; j < mesh.vertices.length; j++)
						for (int k = 0; k < 3; k++)
							mesh.vertexNormals[j][k] = 0.0f;

					// loop over all faces to contribute
					for (int i = 0; i < mesh.triangles.length; i++)
						for (int j = 0; j < 3; j++)
							for (int k = 0; k < 3; k++)
								mesh.vertexNormals[mesh.triangles[i][j]][k] += mesh.triangleNormals[i][k];

					// loop over all vertices to normalize
					for (int j = 0; j < mesh.vertices.length; j++) {
						normalize(mesh.vertexNormals[j]);
					}

					// these are all preparations for phong shading
					if (mesh.fgp == 'p' || mesh.fgp == 'P')
						break;

					// for gouraud-shading: initialize and calculate vertex
					// colors
					mesh.vertexColors = new float[mesh.vertices.length][3];

					rayEx = 0.0f;
					rayEy = 0.0f;
					rayEz = 0.0f;

					// loop over all vertices
					for (int i = 0; i < mesh.vertices.length; i++) {
						// the intersection point is the vertex
						ip = mesh.vertices[i];

						// the material is the material of the vertex
						matIndex = mesh.verticesMat[i];
						material = mesh.materials[matIndex];
						materialN = mesh.materialsN[matIndex];

						// x, y, z: view coordinates are intersection point
						x = ip[0];
						y = ip[1];
						z = ip[2];

						// ray vector
						rayVx = x - rayEx;
						rayVy = y - rayEy;
						rayVz = z - rayEz;

						// fetch precalculated vertex normal
						n = mesh.vertexNormals[i];

						rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

						// backface? => next vertex
						if (rayVn >= 0)
							continue;

						// light vector at the intersection point
						l[0] = ICenter[0] - ip[0];
						l[1] = ICenter[1] - ip[1];
						l[2] = ICenter[2] - ip[2];
						normalize(l);

						// viewing vector at intersection point
						v[0] = -rayVx;
						v[1] = -rayVy;
						v[2] = -rayVz;
						normalize(v);

						// illuminate
						color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

						// write color to vertex
						mesh.vertexColors[i][0] = color.getRed() / 255f;
						mesh.vertexColors[i][1] = color.getGreen() / 255f;
						mesh.vertexColors[i][2] = color.getBlue() / 255f;
					}
				}

			} else if (scene.getHeader() == "Obj") {
				Obj_Mesh obj_mesh = (Obj_Mesh) scene;

				switch (obj_mesh.fgp) {
				case 'f':
				case 'F':
					// for flat-shading: initialize and calculate triangle
					// colors
					obj_mesh.triangleColors = new float[obj_mesh.triangles.length][3];

					rayEx = 0.0f;
					rayEy = 0.0f;
					rayEz = 0.0f;

					// loop over all triangles
					for (int i = 0; i < obj_mesh.triangles.length; i++) {
						// the intersection point is the first vertex of the
						// triangle
						ip = obj_mesh.vertices[obj_mesh.triangles[i][0]];

						// the material is the material of the first triangle
						// point
						matIndex = obj_mesh.verticesMat[obj_mesh.triangles[i][0]];
						material = obj_mesh.materials[matIndex];
						materialN = obj_mesh.materialsN[matIndex];

						// x, y, z: view coordinates are intersection point
						x = ip[0];
						y = ip[1];
						z = ip[2];

						// ray vector
						rayVx = x - rayEx;
						rayVy = y - rayEy;
						rayVz = z - rayEz;

						// fetch precalculated face normal
						n = obj_mesh.triangleNormals[i];

						rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

						// backface? => next triangle
						if (rayVn >= 0)
							continue;

						// light vector at the intersection point
						l[0] = ICenter[0] - ip[0];
						l[1] = ICenter[1] - ip[1];
						l[2] = ICenter[2] - ip[2];
						normalize(l);

						// viewing vector at intersection point
						v[0] = -rayVx;
						v[1] = -rayVy;
						v[2] = -rayVz;
						normalize(v);

						// illuminate
						color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

						// write color to triangle
						obj_mesh.triangleColors[i][0] = color.getRed() / 255f;
						obj_mesh.triangleColors[i][1] = color.getGreen() / 255f;
						obj_mesh.triangleColors[i][2] = color.getBlue() / 255f;

					}

					break;
				}
			} else if (scene.getHeader() == "GLTF") {
				GLTF_Mesh gltf_mesh = (GLTF_Mesh) scene;

				switch (gltf_mesh.fgp) {
				case 'f':
				case 'F':
					// for flat-shading: initialize and calculate triangle
					// colors
					gltf_mesh.triangleColors = new float[gltf_mesh.triangles.length][3];

					rayEx = 0.0f;
					rayEy = 0.0f;
					rayEz = 0.0f;

					// loop over all triangles
					for (int i = 0; i < gltf_mesh.triangles.length; i++) {
						// the intersection point is the first vertex of the
						// triangle
						ip = gltf_mesh.vertices[gltf_mesh.triangles[i][0]];

						// the material is the material of the first triangle
						// point
						matIndex = gltf_mesh.verticesMat[gltf_mesh.triangles[i][0]];
						material = gltf_mesh.materials[matIndex];
						materialN = gltf_mesh.materialsN[matIndex];

						// x, y, z: view coordinates are intersection point
						x = ip[0];
						y = ip[1];
						z = ip[2];

						// ray vector
						rayVx = x - rayEx;
						rayVy = y - rayEy;
						rayVz = z - rayEz;

						// fetch precalculated face normal
						n = gltf_mesh.triangleNormals[i];

						rayVn = rayVx * n[0] + rayVy * n[1] + rayVz * n[2];

						// backface? => next triangle
						if (rayVn >= 0)
							continue;

						// light vector at the intersection point
						l[0] = ICenter[0] - ip[0];
						l[1] = ICenter[1] - ip[1];
						l[2] = ICenter[2] - ip[2];
						normalize(l);

						// viewing vector at intersection point
						v[0] = -rayVx;
						v[1] = -rayVy;
						v[2] = -rayVz;
						normalize(v);

						// illuminate
						color = phongIlluminate(material, materialN, l, n, v, Ia, Ids);

						// write color to triangle
						gltf_mesh.triangleColors[i][0] = color.getRed() / 255f;
						gltf_mesh.triangleColors[i][1] = color.getGreen() / 255f;
						gltf_mesh.triangleColors[i][2] = color.getBlue() / 255f;

					}

					break;
				}
			}
		}

		System.out.println("Vorverarbeitung 2 beendet");
	}

	public static void main(String[] args) {
		Raytracer00 rt = new Raytracer00();

		// rt.doRayTrace();
	}
}
