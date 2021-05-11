using System;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
public class PlayerCollider : MonoBehaviour
{
    public int nrPoints = 18;
    public float colliderOffset = 1.2f;
    public GameObject DistanceSphere;
    public GameObject colPlane;

    private DistanceFunctions Df;
    private PlayerController playCon;
    private float colRadius;
    private RaymarchCamera camScript;
    private Vector3 ro;
    Vector3[] colliders;
    GameObject[] lines;
    GameObject modLine;
    void Start()
    {
        camScript = Camera.main.GetComponent<RaymarchCamera>();
        colRadius = GetComponent<SphereCollider>().radius;
        Df = GetComponent<DistanceFunctions>();
        playCon = GetComponent<PlayerController>();
        colliders = PointsOnSphere(nrPoints, colRadius * colliderOffset);
    }
    void Update()
    {
        if (IsClose())
        {
            RayMarch(colliders);
        }
    }
        // the distancefunction from the player
    public float DistanceField(Vector3 p)
    {
        
        float dist = Camera.main.farClipPlane + 1;

        //menger sponge
        if (camScript._drawMergerSponge)
        {
            dist = Df.sdMerger(p, camScript._GlobalScale, camScript._iterations, camScript._modOffsetPos, camScript._iterationTransform.inverse, camScript._globalTransform.inverse, camScript._smoothRadius, camScript._scaleFactor);

        }
        // menger cylinder
        else if (camScript._drawTowerIFS)
        {
            dist = Df.towerIFS(p, camScript._GlobalScale);
        }

        return dist;
    }
    // the raymarcher from the player
    void RayMarch(Vector3[] rd)
    {
        ro = transform.position;
        int nrHits = 0;

        for (int i = 0; i < rd.Length; i++)
        {
            Vector3 p = ro + Vector3.Normalize(rd[i]) * colRadius;
            //check hit
            float d = DistanceField(p);


            if (d < 0.001) //hit
            {
                Debug.Log("hit" + i);
                nrHits++;
                //collision
                SetColPlane(rd[i]);

            }
            // resets player position if stuck with to many collisions at once
            if (nrHits > rd.Length * 0.45f)
            {
                transform.position = new Vector3(0, 0, camScript._GlobalScale * 1.6f);
            }
            
        }
    }
    // sets the collision plane
    private void SetColPlane(Vector3 hitPoint)
    {
        Instantiate(colPlane, hitPoint + transform.position, Quaternion.identity);
    }
    // checks if the player is close
    bool IsClose()
    {
        float d = DistanceField(transform.position);
        playCon.maxSpeed = Mathf.Min(playCon.maxMaxSpeed, Mathf.Sqrt(d) * playCon.initMaxSpeed); // player speed is regulated by distance
        return d - (colRadius * colliderOffset) < 0.001;
    }
    //creates a fixed number of points on a sphere
    Vector3[] PointsOnSphere(int n, float b)
    {
        List<Vector3> upts = new List<Vector3>();
        float inc = Mathf.PI * (3 - Mathf.Sqrt(5));
        float off = 2.0f / n;
        float x = 0;
        float y = 0;
        float z = 0;
        float r = 0;
        float phi = 0;

        for (var k = 0; k < n; k++)
        {
            y = k * off - 1 + (off / 2);
            r = Mathf.Sqrt(1 - y * y);
            phi = k * inc;
            x = Mathf.Cos(phi) * r;
            z = Mathf.Sin(phi) * r;

            upts.Add(new Vector3(x, y, z) * b);
        }
        Vector3[] pts = upts.ToArray();
        return pts;
    }
}
