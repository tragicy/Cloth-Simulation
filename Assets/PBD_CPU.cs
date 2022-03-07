using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class PBD_CPU : MonoBehaviour
{

	public Transform sphereTrans;
	Vector3 sphereCenter;
	public SphereSimulation sphereSimulation;
	public BoxSimulation boxSimulation;
	public CapsuleSimulation capsuleSimulation;
	float r = 2.5f;
	float t = 0.0333f;
	float damping = 0.99f;
	float Un = 0.5f;
	float Ut = 0.9f;

	int[] E;
	float[] L;
	Vector3[] V;
	Vector3[] X;
	[Range(0,3.0f)]
	public float softness = 0.5f;
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
		//int[] TT = new int[T.Length - 6];
		//for (int i = 0; i < TT.Length; i++)
		//{
		//	TT[i] = T[i];
		//}

		//Vector3[] X = mesh.vertices;
		//Vector2[] UV = mesh.uv;
		//int[] T = mesh.triangles;

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
		List<int> potentialEdge = new List<int>();
		for (int i = 0; i < E.Length/2; i++)
		{
			//int count = 0;
			int ii = E[i * 2 + 0];
			int jj = E[i * 2 + 1];
			List<int> potentialPoint = new List<int>();
			for (int j = 0; j < T.Length/3; j++)
			{
				int containsCount = 0;
				int p0 = T[j * 3 + 0];
				int p1 = T[j * 3 + 1];
				int p2 = T[j * 3 + 2];
				if (ii == p0 || ii == p1 || ii == p2)
					containsCount++;
				if (jj == p0 || jj == p1 || jj == p2)
					containsCount++;
				if (containsCount == 2)
				{
					int p = p0 + p1 + p2 - ii - jj;
					potentialPoint.Add(p);
				}
			}
			if (potentialPoint.Count == 2)
			{
				potentialEdge.Add(potentialPoint[0]);
				potentialEdge.Add(potentialPoint[1]);
			}
		}
		int[] EE = new int[E.Length + potentialEdge.Count];
		for (int i = 0; i < E.Length; i++)
		{
			EE[i] = E[i];
		}
		for (int i = 0; i < potentialEdge.Count; i++)
		{
			EE[i + E.Length] = potentialEdge[i];
		}
        E = new int[EE.Length];
        E = EE;
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
	void ModifyMesh(int edgeIndex)
	{
		int ii = E[2 * edgeIndex];
		int jj = E[2 * edgeIndex + 1];

		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = new Vector3[mesh.vertices.Length + 1];
		for (int i = 0; i < X.Length -1; i++)
		{
			X[i] = mesh.vertices[i];
		}
		X[X.Length - 1] = X[ii];
		int newVertex = X.Length -1;
		int[] TT = mesh.triangles;
		for (int i = 0; i < TT.Length; i+=3)
		{
			int count = 0;
			int a0 = TT[i];
			int a1 = TT[i + 1];
			int a2 = TT[i + 2];
			if (a0 == ii)
				count++;
			if (a0 == jj)
				count++;
			if (a1 == ii)
				count++;
			if (a1 == jj)
				count++;
			if (a2 == ii)
				count++;
			if (a2 == jj)
				count++;
			if (count == 2)
			{
				if (a0 == ii)
					a0 = newVertex;
				if (a1 == ii)
					a1 = newVertex;
				if (a2 == ii)
					a2 = newVertex;
			}
		}
		E[2 * edgeIndex] = newVertex;		
		mesh.triangles = TT;
		mesh.vertices = X;
		mesh.RecalculateNormals();
		Vector3[] VV = new Vector3[V.Length + 1];
		for (int i = 0; i < V.Length; i++)
		{
			VV[i] = V[i];
		}
		VV[VV.Length - 1] = V[ii];
		V = VV;
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
			if (ii >= vertices.Length - 1)
				Debug.Log(ii);
			Vector3 Xij = (vertices[jj] - vertices[ii]);
			float l0 = Xij.magnitude;
			float s = 1.0f / (2.0f + softness / t);
            //if (ii != 0 && ii != 20)
                if (ii != 0 && ii != 20 && ii != XX.Length - 1 && ii != XX.Length - 21)
                    vertices[ii] -= s * (L[i] - l0) * Xij.normalized;

            //if (jj != 0 && jj != 20)
                if (jj != 0 && jj != 20 && jj != XX.Length - 1 && jj != XX.Length - 21)
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
            sum_x[ii] = sum_x[ii] + vertices[ii] - 1.0f/ (2.0f + softness/t ) * (l0 - L[i]) * xij.normalized;
            sum_x[jj] = sum_x[jj] + vertices[jj] + 1.0f / (2.0f + softness / t) * (l0 - L[i]) * xij.normalized;

            sum_n[ii]++;
			sum_n[jj]++;
		}
		for (int i = 0; i < sum_x.Length; i++)
		{
			if (i == 0 || i == 20 || i == sum_x.Length-1 || i== sum_x.Length -21) continue;
			Vector3 newXi = (0.2f * vertices[i] + sum_x[i]) / (0.2f + sum_n[i]);
			vertices[i] = newXi;
		}
	}

	void CollisionWithCapsule(CapsuleSimulation.Capsule capsule)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		List<Vector3> overlapPoints = new List<Vector3>();
		//For every vertex, detect collision and apply impulse if needed.
		//...
		for (int i = 0; i < X.Length; i++)
		{
			if (capsule.intersectionTest(X[i],out Vector3 np))
			{
				float phiX = -Vector3.Magnitude(X[i] - np);
				{
					overlapPoints.Add(X[i]);
					Vector3 N = (np - X[i]).normalized;
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
		}
		mesh.vertices = X;
	}
	void CollisionWithAABB(BoxSimulation.Box box)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		List<Vector3> overlapPoints = new List<Vector3>();
		//For every vertex, detect collision and apply impulse if needed.
		//...
		for (int i = 0; i < X.Length; i++)
		{
			if (box.intersectionTest(X[i]))
			{
				Vector3 np = box.nearPoint(X[i], out Vector3 n);
				float phiX = -Vector3.Magnitude(X[i] - np);
				{
					overlapPoints.Add(X[i]);
					Vector3 N = n;
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
		}
		mesh.vertices = X;
	}
	void CollisionWithSphere()
	{
		sphereCenter = sphereSimulation.ball.x;
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		List<Vector3> overlapPoints = new List<Vector3>();
		//For every vertex, detect collision and apply impulse if needed.
		//...
		for (int i = 0; i < X.Length; i++)
		{
			float distance = Vector3.Magnitude(X[i] - sphereCenter);
			float phiX = distance - r;
			if (phiX < 0.03)
			{
				overlapPoints.Add(X[i]);
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
		for (int i = 0; i < overlapPoints.Count; i++)
		{
			Vector3 N = (sphereCenter - overlapPoints[i]).normalized;
			sphereSimulation.ball.AddForce(1.5f * N);
		}
		mesh.vertices = X;
	}
	void CollisionWithPlan(Vector3 p, Vector3 n)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		for (int i = 0; i < X.Length; i++)
		{
			float phiX =Vector3.Dot((X[i] -p),n);
			if (phiX < 0.03)
			{
				Vector3 N = n;
				X[i] -= N * (phiX - 0.01f);

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
	void SelfIntersectionTest()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		float r = 0.2f;
		for (int i = 0; i < X.Length; i++)
		{
			for (int j = i + 1; j < X.Length; j++)
			{
				Vector3 distance = X[j] - X[i];
				if (distance.magnitude < r + r)
				{
					float corr = r + r - distance.magnitude;
					corr /= 2.0f;

					//update x
					Vector3 dir = distance.normalized;
					X[i] -= dir * corr;
					X[j] += dir * corr;
					float v1 = Vector3.Dot(V[i], dir);
					float v2 = Vector3.Dot(V[j], dir);
					float m1 = 1;
					float m2 = 1;

					float newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * 1.0f) / (m1 + m2);
					float newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * 1.0f) / (m1 + m2);
					V[i] = V[i] + dir * (newV1 - v1);
					V[j] = V[j] + dir * (newV2 - v2);
				}
			}
		}
		//Vector3[] Next_X = new Vector3[X.Length];
		//for (int i = 0; i < X.Length; i++)
		//{
		//	Next_X[i] = X[i] + V[i] * t;
		//}
		//for (int i = 0; i < Next_X.Length; i++)
		//{
		//	for (int j = i + 1; j < Next_X.Length; j++)
		//	{
		//		Vector3 distance = Next_X[j] - Next_X[i];
		//		if (distance.magnitude < r + r)
		//		{
		//			float corr = r + r - distance.magnitude;
		//			corr /= 2.0f;

		//			//update x
		//			Vector3 dir = distance.normalized;
		//			X[i] -= dir * corr;
		//			X[j] += dir * corr;
		//			float v1 = Vector3.Dot(V[i], dir);
		//			float v2 = Vector3.Dot(V[j], dir);
		//			float m1 = 1;
		//			float m2 = 1;

		//			float newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * 1.0f) / (m1 + m2);
		//			float newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * 1.0f) / (m1 + m2);
		//			V[i] = V[i] + dir * (newV1 - v1);
		//			V[j] = V[j] + dir * (newV2 - v2);
		//		}
		//	}
		//}
		mesh.vertices = X;
	}
	// Update is called once per frame
	void Update()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;
		for (int i = 0; i < X.Length; i++)
		{
            if (i == 0 || i == 20 || i == X.Length - 1 || i == X.Length - 21) continue;
            //if (i == 0 || i == 20) continue;

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
        CollisionWithSphere();
        CollisionWithAABB(boxSimulation.box);
        CollisionWithCapsule(capsuleSimulation.capsule);
        //SelfIntersectionTest();
        CollisionWithPlan(new Vector3(0, -5, 0), new Vector3(0, 1, 0));

		mesh.RecalculateNormals();
	}
}
