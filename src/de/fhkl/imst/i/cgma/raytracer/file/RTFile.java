package de.fhkl.imst.i.cgma.raytracer.file;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.Map;

public abstract class RTFile {
	public abstract String getHeader();

	protected abstract void readContent(LineNumberReader f) throws IOException;

	protected String fileName;

	public String getFileName() {
		return fileName;
	}

	@SuppressWarnings("serial")
	public static Map<String, Class<? extends RTFile>> classMapping = new HashMap<String, Class<? extends RTFile>>() {
		{
			put("TRIANGLE_MESH", T_Mesh.class);
			put("IMPLICIT_SPHERE", I_Sphere.class);
			put("SCENE_GRAPH", RTScene.class);
			put("Obj", Obj_Mesh.class);
		}
	};

	public static Class<? extends RTFile> getType(File f) {
		FileReader fr;
		try {
			fr = new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String header = readLine(br);
			System.out.println(header);
			//
			// Dateien des Typs .obj werden in der Klasse Obj_Mesh auf korrektheit getestet
			//
			if (getBaseType(f.getName()).equals("obj"))
				return Obj_Mesh.class;
			if (classMapping.containsKey(header))
				return classMapping.get(header);
			return null;
		} catch (IOException e) {
			System.out.println(e);
			return null;
		}
	}

	/**
	 * Gets the base type, without extension, of given file name.
	 * <p/>
	 * e.g. getBaseType("file.txt") will return "txt"
	 *
	 * @param fileName
	 * @return the file typ
	 */
	public static String getBaseType(String fileName) throws IOException {
		int index = fileName.lastIndexOf('.');
		if (index == -1) {
			throw new IOException("Ungültiger Filetype");
		} else {
			return fileName.substring(index + 1, fileName.length());
		}
	}

	public static <FTYPE extends RTFile> FTYPE read(Class<FTYPE> _class, File f) throws IOException {
		try {
			FTYPE result = _class.newInstance();
			result.fileName = f.getName();
			// check header
			FileReader fr = new FileReader(f);
			LineNumberReader br = new LineNumberReader(fr);
			//
			// Dateien des Typs .obj werden in der Klasse Obj_Mesh auf korrektheit getestet
			//
			if (getBaseType(f.getName()).equals("obj")) {
				// nothing
			} else if (!readLine(br).toLowerCase().equals(result.getHeader().toLowerCase())
					&& !getBaseType(f.getName()).equals("obj"))
				throw new IOException("Ungültiger header");

			result.readContent(br);

			return result;
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return null;
	}

	protected static String readLine(BufferedReader br) throws IOException {
		String result = br.readLine().trim();
		while (result.startsWith("#") || result.trim().isEmpty())
			result = br.readLine().trim();
		return result;
	}
}
