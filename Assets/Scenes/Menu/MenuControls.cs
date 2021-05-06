using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
public class MenuControls : MonoBehaviour
{
    // Start is called before the first frame update
    public void SetInterier()
    {
        SceneManager.LoadScene("Interier");
    }
    public void SetTown()
    {
        SceneManager.LoadScene("LSystemTown");
    }
    public void SetLandscape()
    {
        SceneManager.LoadScene("Landscape");
    }
    public void ExitPressed()
    {
        Debug.Log("Exit pressed!");
        Application.Quit();
    }
}
