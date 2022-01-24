using UnityEngine;
using System.Collections;

public class PBD_model: MonoBehaviour {

	public ComputeShader shader;
	public Transform sphereTrans;
	float r = 2.5f;
	Vector3 sphereCenter;
	float 		t= 0.0333f;
	float		damping= 0.99f;
	float Un = 0.5f;
	float Ut = 0.9f;
	static int n = 21;
	int SIZE = n * n;

	int[] 		E;
	float[] 	L;
	int[] sum_n;

	Vector3[] V;
	Vector3[] Sum_X;
	Vector3[] X;
	//Vector3[]	D;
	Vector3 g = new Vector3(0, -9.80f, 0);


	#region ComputeShader
	int kernelHandle;
	int kernelHandle1;

	ComputeBuffer XBuffer;
	ComputeBuffer Sum_xBuffer;
	ComputeBuffer VBuffer;
	ComputeBuffer EBuffer;
	ComputeBuffer LBuffer;
	ComputeBuffer Sum_nBuffer;
	void InitData()
	{
		kernelHandle = shader.FindKernel("CSMain");
		kernelHandle1 = shader.FindKernel("CSMain1");

		XBuffer = new ComputeBuffer(SIZE, 3 * sizeof(float));
		XBuffer.SetData(X);
		shader.SetBuffer(kernelHandle, "X", XBuffer);
		shader.SetBuffer(kernelHandle1, "X", XBuffer);

		Sum_xBuffer = new ComputeBuffer(SIZE, 3 * sizeof(float));
		Sum_xBuffer.SetData(Sum_X);
		shader.SetBuffer(kernelHandle, "Sum_X", Sum_xBuffer);
		shader.SetBuffer(kernelHandle1, "Sum_X", Sum_xBuffer);

		VBuffer = new ComputeBuffer(SIZE, 3 * sizeof(float));
		VBuffer.SetData(V);
		shader.SetBuffer(kernelHandle, "V", VBuffer);
		shader.SetBuffer(kernelHandle1, "V", VBuffer);

		EBuffer = new ComputeBuffer(E.Length, sizeof(int));
		EBuffer.SetData(E);
		shader.SetBuffer(kernelHandle, "E", EBuffer);
		shader.SetBuffer(kernelHandle1, "E", EBuffer);

		LBuffer = new ComputeBuffer(L.Length, sizeof(float));
		LBuffer.SetData(L);
		shader.SetBuffer(kernelHandle, "L", LBuffer);
		shader.SetBuffer(kernelHandle1, "L", LBuffer);

		Sum_nBuffer = new ComputeBuffer(sum_n.Length, sizeof(int));
		Sum_nBuffer.SetData(sum_n);
		shader.SetBuffer(kernelHandle, "Sum_n", Sum_nBuffer);
		shader.SetBuffer(kernelHandle1, "Sum_n", Sum_nBuffer);
	}
	#endregion


	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		//Resize the mesh.
		X  	= new Vector3[n*n];
		Vector2[] UV 	= new Vector2[n*n];
		int[] T	= new int[(n-1)*(n-1)*6];
		for(int j=0; j<n; j++)
		for(int i=0; i<n; i++)
		{
			X[j*n+i] =new Vector3(5-10.0f*i/(n-1), 0, 5-10.0f*j/(n-1));
			UV[j*n+i]=new Vector3(i/(n-1.0f), j/(n-1.0f));
		}
		int t=0;
		for(int j=0; j<n-1; j++)
		for(int i=0; i<n-1; i++)	
		{
			T[t*6+0]=j*n+i;
			T[t*6+1]=j*n+i+1;
			T[t*6+2]=(j+1)*n+i+1;
			T[t*6+3]=j*n+i;
			T[t*6+4]=(j+1)*n+i+1;
			T[t*6+5]=(j+1)*n+i;
			t++;
		}
		mesh.vertices	= X;
		mesh.triangles	= T;
		mesh.uv 		= UV;
		mesh.RecalculateNormals ();

		//Construct the original edge list
		int[] _E = new int[T.Length*2];
		for (int i=0; i<T.Length; i+=3) 
		{
			_E[i*2+0]=T[i+0];
			_E[i*2+1]=T[i+1];
			_E[i*2+2]=T[i+1];
			_E[i*2+3]=T[i+2];
			_E[i*2+4]=T[i+2];
			_E[i*2+5]=T[i+0];
		}
		//Reorder the original edge list
		for (int i=0; i<_E.Length; i+=2)
			if(_E[i] > _E[i + 1]) 
				Swap(ref _E[i], ref _E[i+1]);
		//Sort the original edge list using quicksort
		Quick_Sort (ref _E, 0, _E.Length/2-1);

		int e_number = 0;
		for (int i=0; i<_E.Length; i+=2)
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
				e_number++;

		E = new int[e_number * 2];
		for (int i=0, e=0; i<_E.Length; i+=2)
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
			{
				E[e*2+0]=_E [i + 0];
				E[e*2+1]=_E [i + 1];
				e++;
			}

		L = new float[E.Length/2];
		//D = new Vector3[E.Length / 2];
		for (int e=0; e<E.Length/2; e++) 
		{
			int i = E[e*2+0];
			int j = E[e*2+1];
			L[e]=(X[i]-X[j]).magnitude;
		}

		V = new Vector3[X.Length];
		for (int i=0; i<X.Length; i++)
			V[i] = new Vector3 (0, 0, 0);

		sum_n = new int[mesh.vertices.Length];
		Sum_X = new Vector3[mesh.vertices.Length];
		Debug.Log(X.Length);
		//Debug.Log(E.Length);
		//Debug.Log(L.Length);

		InitData();
    }

	void Quick_Sort(ref int[] a, int l, int r)
	{
		int j;
		if(l<r)
		{
			j=Quick_Sort_Partition(ref a, l, r);
			Quick_Sort (ref a, l, j-1);
			Quick_Sort (ref a, j+1, r);
		}
	}

	int  Quick_Sort_Partition(ref int[] a, int l, int r)
	{
		int pivot_0, pivot_1, i, j;
		pivot_0 = a [l * 2 + 0];
		pivot_1 = a [l * 2 + 1];
		i = l;
		j = r + 1;
		while (true) 
		{
			do ++i; while( i<=r && (a[i*2]<pivot_0 || a[i*2]==pivot_0 && a[i*2+1]<=pivot_1));
			do --j; while(  a[j*2]>pivot_0 || a[j*2]==pivot_0 && a[j*2+1]> pivot_1);
			if(i>=j)	break;
			Swap(ref a[i*2], ref a[j*2]);
			Swap(ref a[i*2+1], ref a[j*2+1]);
		}
		Swap (ref a [l * 2 + 0], ref a [j * 2 + 0]);
		Swap (ref a [l * 2 + 1], ref a [j * 2 + 1]);
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

		//Apply PBD here.
		//...
		//for (int i = 0; i < L.Length; i++)
		//{
		//    int ii = E[2 * i];
		//    int jj = E[2 * i + 1];
		//    Vector3 Xij = (X[ii] - X[jj]).normalized;
		//    Sum_X[ii] = Sum_X[ii] + 0.5f * (X[ii] + X[jj] + L[i] * Xij);
		//    Sum_X[jj] = Sum_X[jj] + 0.5f * (X[ii] + X[jj] - L[i] * Xij);
		//    sum_n[ii]++;
		//    sum_n[jj]++;
		//}
		for (int k = 0; k < 32; k++)
		{
			for (int i = 0; i < Sum_X.Length; i++)
			{
				for (int j = 0; j < L.Length; j++)
				{
					if (E[2 * j] == i)
					{
						int ii = E[2 * j];
						int jj = E[2 * j + 1];
						Vector3 Xij = (X[ii] - X[jj]).normalized;
						Sum_X[ii] = Sum_X[ii] + 0.5f * (X[ii] + X[jj] + L[j] * Xij);
						sum_n[ii]++;
					}
					else if (E[2 * j + 1] == i)
					{
						int ii = E[2 * j];
						int jj = E[2 * j + 1];
						Vector3 Xij = (X[ii] - X[jj]).normalized;
						Sum_X[jj] = Sum_X[jj] + 0.5f * (X[ii] + X[jj] - L[j] * Xij);
						sum_n[jj]++;
					}
				}
				if (i == 0 || i == 20) continue;

				Vector3 newXi = (0.5f * X[i] + Sum_X[i]) / (0.5f + sum_n[i]);
				V[i] = V[i] + (newXi - X[i]) / t;
				X[i] = newXi;
				Sum_X[i] = Vector3.zero;
				sum_n[i] = 0;
			}
		}
    }

	void Collision_Handling()
	{
		sphereCenter = sphereTrans.position;
		//For every vertex, detect collision and apply impulse if needed.
		//...
		for (int i = 0; i < X.Length; i++)
		{
			float distance = Vector3.Magnitude(X[i] - sphereCenter);
			float phiX = distance - r;
			if (phiX < 0.03)
			{
				Vector3 N = (X[i] - sphereCenter).normalized;
				X[i] -= N * (phiX - 0.05f);

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

	// Update is called once per frame
	int num = 0;
	void Update () 
	{
		num++;
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		X = mesh.vertices;
		for(int i=0; i<X.Length; i++)
		{
			if(i<=20)	continue;
			//Initial Setup
			//...
			V[i] = V[i] + g * t;
			X[i] = X[i] + V[i] * t;
			V[i] *= damping;
		}

		//for(int l=0; l<32; l++)
			//Strain_Limiting ();

		
		XBuffer.SetData(X);
		VBuffer.SetData(V);
		//Sum_xBuffer.SetData(Sum_X);
		//Sum_nBuffer.SetData(sum_n);
		shader.SetBuffer(kernelHandle1, "X", XBuffer);
		shader.SetBuffer(kernelHandle1, "V", VBuffer);
		//shader.SetBuffer(kernelHandle1, "Sum_X", Sum_xBuffer);
		//shader.SetBuffer(kernelHandle1, "Sum_n", Sum_nBuffer);
		
		//shader.Dispatch(kernelHandle, L.Length, 1, 1);
		for(int i =0; i<16;i++)
		shader.Dispatch(kernelHandle1, SIZE, 1, 1);
		XBuffer.GetData(X);
		VBuffer.GetData(V);
		Collision_Handling();

		mesh.vertices = X;




		mesh.RecalculateNormals();
		//sum_n = new int[mesh.vertices.Length];
		//Sum_X = new Vector3[mesh.vertices.Length];



		//XBuffer.GetData(X);
		//VBuffer.GetData(V);

		//V = V;
		//mesh.vertices = X;
		//mesh.RecalculateNormals();

	}


}

