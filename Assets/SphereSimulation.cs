using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SphereSimulation : MonoBehaviour
{
    public GameObject SpherePrefab;
    public class Ball
    {
        GameObject go;
        public Vector3 v;
        public Vector3 x;
        float dt = 1.0f / 60.0f;
        public float mass = 2.0f;
        public float radius = 2.5f;
        Vector3 g = new Vector3(0, -10.0f, 0);
        Vector3 Force = Vector3.zero;
        float Un = 0.5f;
        float Ut = 0.9f;
        public Ball(Vector3 _x, Vector3 _v, GameObject goPrefab)
        {
            v = _v;
            x = _x;
            go = Instantiate(goPrefab);
            go.transform.position = x;
            go.SetActive(true);
        }
        public void AddForce(Vector3 force)
        {
            Force += force;
        }
        public void Simulate()
        {
            AddForce(g * mass);
            Vector3 a = Force / mass;
            v += a * dt;
            x += v * dt;
            go.transform.position = x;
            Force = Vector3.zero;
        }
        public void CollisionWithFloor(Vector3 p, Vector3 n)
        {
            float phiX = Vector3.Dot((x - p), n);
            if (phiX < radius)
            {
                //overlapPoints.Add(X[i]);
                Vector3 N = n;
                x -= N * (phiX - 0.01f - radius);

                Vector3 _v = v;
                if (Vector3.Dot(_v, N) < 0.0f)
                {
                    Vector3 Vn = Vector3.Dot(_v, N) * N;
                    Vector3 Vt = _v - Vn;
                    Vector3 VNewN = -Un * Vn;

                    float a = 1 - Ut * (1 + Un) * Vn.magnitude / Vt.magnitude;
                    if (a < 0)
                        a = 0;
                    Vector3 VNewt = a * Vt;
                    Vector3 VNew = VNewt + VNewN;
                    v = VNew;
                }
            }
        }
    }
    public Ball ball;
    // Start is called before the first frame update
    void Start()
    {
        ball = new Ball(new Vector3(0,-1,0), Vector3.zero, SpherePrefab);
    }

    // Update is called once per frame
    void Update()
    {
        ball.Simulate();
        ball.CollisionWithFloor(new Vector3(0, -5, 0), new Vector3(0, 1, 0));
    }
}
