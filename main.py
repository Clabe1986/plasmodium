import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle
import streamlit as st
import subprocess
import csv
from pypdb import *


@st.cache_data
def calculate_lipinski_descriptors(smiles):
    """
    Calculate Lipinski descriptors for a given SMILES string.

    Parameters:
        smiles (str): Canonical SMILES string of a compound.

    Returns:
        pandas.DataFrame: DataFrame containing the calculated Lipinski descriptors.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    descriptors = {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol)
    }

    return pd.DataFrame(descriptors, index=[0])


def generate_csv_file(string1, string2, filename):
    """
    Generate a CSV file with a single row containing the concatenated string of string1 and string2.

    Parameters:
        string1 (str): First string.
        string2 (str): Second string.
        filename (str): Name of the output CSV file.
    """
    data = [[string1 + '\t' + string2]]  # Create a list of lists containing the strings

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)


# Set page title and initial layout
st.set_page_config(page_title="Plasmodium Drug Prediction System", page_icon='🦟', layout="wide")

# Apply custom CSS to modify app appearance
st.markdown(
    """
    <style>
    .title {
        color: #FF0000;
        font-size: 32px;
        margin-bottom: 20px;
    }

    .result-text {
        font-size: 20px;
        margin-top: 10px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# App title
st.title("A Plasmodium faliparum Drug Prediction System")

# ---Use local CSS---
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


local_css("style/style.css")


# Main function
def main():
    """
    Main function for the Plasmodium faliparum Drug Prediction System web application.
    """
    # Text input
    smiles_input = st.text_input("Enter Canonical SMILES")

    # Selectbox options
    options = [
        "Compute Lipinski's Descriptors",
        "Predict the Compound's Activity",
        "Predict the Compound's pIC50",
        "Retrieve interacting proteins"
    ]

    # Selectbox
    selected_option = st.selectbox("Select the program to run", options)

    # Button to run the selected option
    if st.button("Run"):
        if not smiles_input:
            st.text("Enter Canonical SMILES")
        else:
            try:
                if selected_option == "Compute Lipinski's Descriptors":
                    # Calculate Lipinski's descriptors and display the result
                    df = calculate_lipinski_descriptors(smiles_input)
                    df.columns = ['Molecular Weight', 'Octanol-Water Partition Coefficient (LogP)',
                                  'Number of Hydrogen Bond Donors', 'Number of Hydrogen Bond Acceptors']
                    hide_table_row_index = """
                        <style>
                        thead tr th:first-child {display:none}
                        tbody th {display:none}
                        </style>
                    """
                    st.markdown(hide_table_row_index, unsafe_allow_html=True)
                    st.table(df)
                elif selected_option == "Predict the Compound's Activity":
                    # Load the model and predict the compound's activity
                    filename = 'lipinsky_model.pkl'
                    loaded_model = pickle.load(open(filename, 'rb'))
                    y_pred = loaded_model.predict(calculate_lipinski_descriptors(smiles_input))
                    if y_pred == [1]:
                        st.text('Active')
                    else:
                        st.text('Inactive')
                elif selected_option == "Predict the Compound's pIC50":
                    # Generate a CSV file and compute pIC50
                    string1 = smiles_input
                    string2 = 'Compound_name'

                    filename = "molecule.smi"
                    generate_csv_file(string1, string2, filename)

                    script_path = "./padel.sh"
                    subprocess.call(['bash', script_path])

                    data = pd.read_csv('descriptors_output.csv')
                    X = data.drop(columns=['Name'])

                    loaded_model = pickle.load(open('model.pkl', 'rb'))
                    y_pred = loaded_model.predict(X)
                    predicted_value = y_pred[0]
                    predicted_value = format(predicted_value, ".2f")
                    st.text("The pIC50 of your compound is " + str(predicted_value))
                elif selected_option == "Retrieve interacting proteins":
                    # Retrieve interacting proteins from PDB database
                    result = Query(smiles_input).search()
                    ids = []
                    for pdb_id in result:
                        if len(ids) < 10:
                            ids.append(pdb_id)
                    ids_str = ", ".join(ids)
                    st.markdown(ids_str)

            except ValueError as e:
                st.error(str(e))
    # ---- CONTACT ----
    with st.container():
        st.write("---")
        st.header("Please contact us!")

        # Documention: https://formsubmit.co/ !!! CHANGE EMAIL ADDRESS !!!
        contact_form = """
          <form action="https://formsubmit.co/Wacheed21@gmail.com" method="POST">
              <input type="hidden" name="_captcha" value="false">
              <input type="text" name="name" placeholder="Your name" required>
              <input type="email" name="email" placeholder="Your email" required>
              <textarea name="message" placeholder="Your message here" required></textarea>
              <button type="submit">Send</button>
          </form>
          """
        left_column, right_column = st.columns(2)
        with left_column:
            st.markdown(contact_form, unsafe_allow_html=True)
        with right_column:
            st.empty()


if __name__ == "__main__":
    main()
