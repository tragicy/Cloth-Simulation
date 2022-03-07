using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BoxSimulation : MonoBehaviour
{
    public class Box
    {
       public Vector3 center;
       public float x_len;
       public float y_len;
       public float z_len;

        public Box(Vector3 _c,float _x, float _y, float _z)
        {
            center = _c;
            x_len = _x;
            y_len = _y;
            z_len = _z;
        }
        public bool intersectionTest(Vector3 point)
        {
            if (point.x < center.x - x_len || point.x > center.x + x_len)
                return false;
            if (point.y < center.y - y_len || point.y > center.y + y_len)
                return false;
            if (point.z < center.z - z_len || point.z > center.z + z_len)
                return false;
            return true;
        }

        public Vector3 nearPoint(Vector3 point,out Vector3 n)
        {
            Vector3 np = Vector3.zero;
            n = Vector3.zero;
            np.x = Mathf.Abs(center.x + x_len - point.x) < Mathf.Abs(center.x - x_len - point.x) ? center.x + x_len - point.x : (center.x - x_len) - point.x;
            np.y = Mathf.Abs(center.y + y_len - point.y) < Mathf.Abs(center.y - y_len - point.y) ? center.y + y_len - point.y : (center.y - y_len) - point.y;
            np.z = Mathf.Abs(center.z + z_len - point.z) < Mathf.Abs(center.z - z_len - point.z) ? center.z + z_len - point.z : (center.z - z_len) - point.z;
            Vector3 tmpNp = np;
            if (Mathf.Abs(np.x) < Mathf.Abs(np.y))
            {
                np.y = 0;
                if (Mathf.Abs(np.x) < Mathf.Abs(np.z))
                    np.z = 0;
                else
                    np.x = 0;
            }
            else
            {
                np.x = 0;
                if (Mathf.Abs(np.y) <= Mathf.Abs(np.z))
                    np.z = 0;
                else
                    np.y = 0;
            }
            float distance = np.magnitude;
            float buffer = 1.0f;
            if (Mathf.Abs(distance - tmpNp.x) < buffer)
            {
                np.x = tmpNp.x;
            }
            if (Mathf.Abs(distance - tmpNp.y) < buffer)
            {
                np.y = tmpNp.y;
            }
            if (Mathf.Abs(distance - tmpNp.z) < buffer)
            {
                np.z = tmpNp.z;
            }
            n = np.normalized;

            np += point;
            return np;
        }
    }
    public Box box;
    // Start is called before the first frame update
    void Start()
    {
       box = new Box(Vector3.zero,2,2,2);
    }

    // Update is called once per frame
    void Update()
    {
        box.center = transform.position;
    }
}
