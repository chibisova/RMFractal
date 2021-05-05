using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;
public class IFSMenu : MonoBehaviour
{
    // Start is called before the first frame update
    public void PlayPressed()
    {
        SceneManager.LoadScene("Menu");
    }
    public void ExitPressed()
    {
        Debug.Log("Exit pressed!");
        Application.Quit();
    }
}