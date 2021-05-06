using UnityEngine;
using System.Collections;
using UnityEngine.UI;

public class RecursiveInstantiator : MonoBehaviour {

	public Canvas canvasFields;
    public int recurse = 5;
    public int splitNumber = 2;
    public Vector3 pivotPosition;

    private InputField[] uiFields;
	// Use this for initialization
	void Awake()
	{
		if (canvasFields != null)
		{
			uiFields = canvasFields.GetComponentsInChildren<InputField>();
		}

	}
	void Start () {
		PopulateFields();
		UpdateParams();

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
	
	// Update is called once per frame
	void Update () {
	
	}
	

	private void UpdateParams()
	{
		if (uiFields != null && uiFields.Length > 0)
		{
			foreach (var field in uiFields)
			{
				if (field.text.Length > 0 && field.name.ToLower().StartsWith("iterations"))
				{
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
	}
	private void PopulateFields()
	{
		if (uiFields != null && uiFields.Length > 0)
		{
			foreach (var field in uiFields)
			{
				if (field.name.ToLower().StartsWith("iterations"))
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
