using UnityEngine;
using System.Collections;
using UnityEngine.UI;

public class RecursiveInstantiator : MonoBehaviour {
	public Canvas canvasFields;

    public int recurse = 8;
    public int splitNumber = 2;
    public Vector3 pivotPosition;
    private InputField[] uiFields;
    void Awake()
    {
	    canvasFields = GameObject.FindGameObjectWithTag("MainCanvas").GetComponent<Canvas>();;
	    Debug.Log("Got" );

	    if (canvasFields != null)
	    {
		    Debug.Log("Inside" );

		    uiFields = canvasFields.GetComponentsInChildren<InputField>();
	    }

    }

	// Use this for initialization
	void Start ()
	{
		PopulateFields();
		Create();
	}

	public void Create()
	{
		recurse -= 1;

		for (int i = 0; i < splitNumber; ++i)
		{
			if (recurse > 0)
			{
				var copy = Instantiate(gameObject);
				var recursive = copy.GetComponent<RecursiveInstantiator>();
				recursive.SendMessage("Generated", new RecursiveBundle() { Index = i, Parent = this });
				
			}
		}
	}
	// ------- UI Fields params -------
	void Update()
	{
	}
	public void TaskOnClick()
	{
		Debug.Log("Set" );
		if (uiFields != null && uiFields.Length > 0)
		{
			Debug.Log("Inside" );
			foreach (var field in uiFields)
			{
				if (field.text.Length > 0 && field.name.ToLower().StartsWith("iter"))
				{
					Debug.Log("Iter" );
					var value = int.Parse(field.text);
					if (recurse != value)
					{
						recurse = value;
					}
				}

				if (field.text.Length > 0 && field.name.ToLower().StartsWith("split"))
				{
					var value = int.Parse(field.text);

					if (splitNumber != value)
					{
						splitNumber = value;
					}
				}
			}
		}
		
		Debug.Log("Pressed!");
		Create();
		
	}


	private void PopulateFields()
	{
		if (uiFields != null && uiFields.Length > 0)
		{
			foreach (var field in uiFields)
			{
				if (field.name.ToLower().StartsWith("iter"))
				{
					field.text = recurse.ToString();
				}

				if (field.name.ToLower().StartsWith("split"))
				{
					field.text = splitNumber.ToString();
				}
			}
		}
	}
	
}
