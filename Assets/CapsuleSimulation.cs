using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CapsuleSimulation : MonoBehaviour
{
    public class Capsule
    {
        public Vector3 pLow;
        public Vector3 pUp;
        float radius;
        public Capsule(Vector3 _pLow, Vector3 _pUp, float r)
        {
            pLow = _pLow;
            pUp = _pUp;
            radius = r;
        }
        public bool projectionTest(Vector3 p)
        {
            Vector3 p01 = pUp - pLow;
            Vector3 p02 = p - pLow;
            Vector3 v = p01.normalized;
            float dotProduct = Vector3.Dot(p02, v);
            if (dotProduct <= p01.magnitude && dotProduct>0)
                return true;
            return false;
        }
        public bool intersectionTest(Vector3 p,out Vector3 np)
        {
            Vector3 p01 = pUp - pLow;
            Vector3 p02 = p - pLow;
            Vector3 p12 = p - pUp;
            Vector3 v = p01.normalized;
            np = p;
            float dotProduct = Vector3.Dot(p02, v);
            //Test lower sphere
            if (dotProduct < 0)
            {
                float distance = p02.magnitude;
                if (distance < radius)
                {
                    np = pLow + p02.normalized * radius;
                    return true;
                }
                else
                    return false;
            }
            else
            {
                //Test middle
                if (dotProduct <= p01.magnitude)
                {
                    Vector3 s = pLow + dotProduct * v;
                    float distance = (p - s).magnitude;
                    if (distance < radius)
                    {
                        np = s + (p - s).normalized * radius;
                        return true;
                    }
                    else
                        return false;
                }
                //Test upper sphere
                else
                {
                    float distance = p12.magnitude;
                    if (distance < radius)
                    {
                        np = pUp + p12.normalized * radius;
                        return true;
                    }
                    else
                        return false;
                }
            }
            //return false;
        }
    }
    public Capsule capsule;
    // Start is called before the first frame update
    void Start()
    {
        capsule = new Capsule(Vector3.zero, new Vector3(0, 2, 0), 2);
        Debug.Log(capsule.projectionTest(new Vector3(1,2.1f,0)));
    }

    // Update is called once per frame
    void Update()
    {
        capsule.pUp = transform.position + new Vector3(0, 2, 0);
        capsule.pLow = transform.position + new Vector3(0, -2, 0);
    }
}
