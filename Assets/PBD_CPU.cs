using UnityEngine;
using System.Collections;

public class PBD_CPU : MonoBehaviour
{

	public Transform sphereTrans;
	float r = 2.5f;
	Vector3 sphereCenter;
	float t = 0.0333f;
	float damping = 0.99f;
	float Un = 0.5f;
	float Ut = 0.9f;

	int[] E;
	float[] L;
	Vector3[] V;
	Vector3[] X;
	float lambda = 0.1f;

	Vector3 g = new Vector3(0, -9.80f, 0);



	// Use this for initialization
	void Start()
	{
		Application.targetFrameRate = 120;
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		//Resize the mesh.
		int n = 21;
		Vector3[] X = new Vector3[n * n];
		Vector2[] UV = new Vector2[n * n];
		int[] T = new int[(n - 1) * (n - 1) * 6];
		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++)
			{
				X[j * n + i] = new Vector3(5 - 10.0f * i / (n - 1), 0, 5 - 10.0f * j / (n - 1));
				UV[j * n + i] = new Vector3(i / (n - 1.0f), j / (n - 1.0f));
			}
		int t = 0;
		for (int j = 0; j < n - 1; j++)
			for (int i = 0; i < n - 1; i++)
			{
				T[t * 6 + 0] = j * n + i;
				T[t * 6 + 1] = j * n + i + 1;
				T[t * 6 + 2] = (j + 1) * n + i + 1;
				T[t * 6 + 3] = j * n + i;
				T[t * 6 + 4] = (j + 1) * n + i + 1;
				T[t * 6 + 5] = (j + 1) * n + i;
				t++;
			}
		mesh.vertices = X;
		mesh.triangles = T;
		mesh.uv = UV;
		mesh.RecalculateNormals();

		//Construct the original edge list
		int[] _E = new int[T.Length * 2];
		for (int i = 0; i < T.Length; i += 3)
		{
			_E[i * 2 + 0] = T[i + 0];
			_E[i * 2 + 1] = T[i + 1];
			_E[i * 2 + 2] = T[i + 1];
			_E[i * 2 + 3] = T[i + 2];
			_E[i * 2 + 4] = T[i + 2];
			_E[i * 2 + 5] = T[i + 0];
		}
		//Reorder the original edge list
		for (int i = 0; i < _E.Length; i += 2)
			if (_E[i] > _E[i + 1])
				Swap(ref _E[i], ref _E[i + 1]);
		//Sort the original edge list using quicksort
		Quick_Sort(ref _E, 0, _E.Length / 2 - 1);

		int e_number = 0;
		for (int i = 0; i < _E.Length; i += 2)
			if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
				e_number++;

		E = new int[e_number * 2];
		for (int i = 0, e = 0; i < _E.Length; i += 2)
			if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1])
			{
				E[e * 2 + 0] = _E[i + 0];
				E[e * 2 + 1] = _E[i + 1];
				e++;
			}

		L = new float[E.Length / 2];
		//D = new Vector3[E.Length / 2];
		for (int e = 0; e < E.Length / 2; e++)
		{
			int i = E[e * 2 + 0];
			int j = E[e * 2 + 1];
			L[e] = (X[i] - X[j]).magnitude;
		}

		V = new Vector3[X.Length];
		for (int i = 0; i < X.Length; i++)
			V[i] = new Vector3(0, 0, 0);

		X = mesh.vertices;
		Debug.Log(E.Length);
		Debug.Log(L.Length);
	}

	void Quick_Sort(ref int[] a, int l, int r)
	{
		int j;
		if (l < r)
		{
			j = Quick_Sort_Partition(ref a, l, r);
			Quick_Sort(ref a, l, j - 1);
			Quick_Sort(ref a, j + 1, r);
		}
	}

	int Quick_Sort_Partition(ref int[] a, int l, int r)
	{
		int pivot_0, pivot_1, i, j;
		pivot_0 = a[l * 2 + 0];
		pivot_1 = a[l * 2 + 1];
		i = l;
		j = r + 1;
		while (true)
		{
			do ++i; while (i <= r && (a[i * 2] < pivot_0 || a[i * 2] == pivot_0 && a[i * 2 + 1] <= pivot_1));
			do --j; while (a[j * 2] > pivot_0 || a[j * 2] == pivot_0 && a[j * 2 + 1] > pivot_1);
			if (i >= j) break;
			Swap(ref a[i * 2], ref a[j * 2]);
			Swap(ref a[i * 2 + 1], ref a[j * 2 + 1]);
		}
		Swap(ref a[l * 2 + 0], ref a[j * 2 + 0]);
		Swap(ref a[l * 2 + 1], ref a[j * 2 + 1]);
		return j;
	}

	void Swap(ref int a, ref int b)
	{
		int temp = a;
		a = b;
		b = temp;
	}

	void Strain_Limiting()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		//Apply PBD here.
		//...
		Vector3[] X = new Vector3[vertices.Length];
		int[] sum_n = new int[vertices.Length];
		for (int i = 0; i < L.Length; i++)
		{
			int ii = E[2 * i];
			int jj = E[2 * i + 1];
			Vector3 Xij = (vertices[ii] - vertices[jj]).normalized;
			X[ii] = X[ii] + 0.5f * (vertices[ii] + vertices[jj] + L[i] * Xij);
			X[jj] = X[jj] + 0.5f * (vertices[ii] + vertices[jj] - L[i] * Xij);
			sum_n[ii]++;
			sum_n[jj]++;
		}
		for (int i = 0; i < X.Length; i++)
		{
			if (i == 0 || i == 20) continue;

			Vector3 newXi = (0.2f * vertices[i] + X[i]) / (0.2f + sum_n[i]);
			V[i] = V[i] + (newXi - vertices[i]) / t;
			vertices[i] = newXi;
		}
		mesh.vertices = vertices;
	}
	void GS_ConstraintProjection(Vector3[] XX)
	{
		Vector3[] vertices = XX;
		//Apply PBD here.
		for (int i = 0; i < L.Length; i++)
		{
			int ii = E[2 * i];
			int jj = E[2 * i + 1];

			Vector3 Xij = (vertices[jj] - vertices[ii]);
			float l0 = Xij.magnitude;
			float s = 1.0f / (2.0f + lambda / t);
			if (ii != 0 && ii != 20)
				vertices[ii] -= s * (L[i] - l0) * Xij.normalized;

			if (jj != 0 && jj != 20)
				vertices[jj] += s * (L[i] - l0) * Xij.normalized;
		}
	}
	void Jacobi_ConstraintProjection(Vector3[] XX)
	{
		Vector3[] vertices = XX;
		//Apply PBD here.
		//...
		//float lambda = 0.0f;
		Vector3[] sum_x = new Vector3[vertices.Length];
		int[] sum_n = new int[vertices.Length];
		for (int i = 0; i < L.Length; i++)
		{
			int ii = E[2 * i];
			int jj = E[2 * i + 1];
            Vector3 xij = vertices[ii] - vertices[jj];
            float l0 = xij.magnitude;
            sum_x[ii] = sum_x[ii] + vertices[ii] - 1.0f/ (2.0f + lambda/t ) * (l0 - L[i]) * xij.normalized;
            sum_x[jj] = sum_x[jj] + vertices[jj] + 1.0f / (2.0f + lambda / t) * (l0 - L[i]) * xij.normalized;

            sum_n[ii]++;
			sum_n[jj]++;
		}
		for (int i = 0; i < sum_x.Length; i++)
		{
			if (i == 0 || i == 20) continue;
			Vector3 newXi = (0.2f * vertices[i] + sum_x[i]) / (0.2f + sum_n[i]);
			vertices[i] = newXi;
		}
	}
	void Collision_Handling()
	{
		sphereCenter = sphereTrans.position;
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		//For every vertex, detect collision and apply impulse if needed.
		//...
		for (int i = 0; i < X.Length; i++)
		{
			float distance = Vector3.Magnitude(X[i] - sphereCenter);
			float phiX = distance - r;
			if (phiX < 0.03)
			{
				Vector3 N = (X[i] - sphereCenter).normalized;
				X[i] -= N * (phiX - 0.1f);

				Vector3 v = V[i];
				if (Vector3.Dot(v, N) < 0)
				{
					Vector3 Vn = Vector3.Dot(v, N) * N;
					Vector3 Vt = v - Vn;
					Vector3 VNewN = -Un * Vn;

					float a = 1 - Ut * (1 + Un) * Vn.magnitude / Vt.magnitude;
					if (a < 0)
						a = 0;
					Vector3 VNewt = a * Vt;
					Vector3 VNew = VNewt + VNewN;
					V[i] = VNew;
				}
			}
		}
		mesh.vertices = X;
	}

	// Update is called once per frame
	void Update()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		for (int i = 0; i < X.Length; i++)
		{
			if (i == 0 || i == 20) continue;
			//Initial Setup
			//...
			V[i] = V[i] + g * t;
			X[i] = X[i] + V[i] * t;
			V[i] *= damping;
		}

		for (int l = 0; l < 16; l++)
		{
			//Jacobi_ConstraintProjection(X);
			GS_ConstraintProjection(X);
		}
		Vector3[] pre_X = mesh.vertices;

		for (int i = 0; i < X.Length; i++)
		{
			V[i] = (X[i] - pre_X[i]) / t;
		}
		mesh.vertices = X;
		Collision_Handling();
		mesh.RecalculateNormals();

	}


}
