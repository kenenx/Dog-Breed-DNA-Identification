from flask import Flask, render_template, request, redirect, url_for, jsonify
from utils import load_sequences
from alignment import find_closest_match
import threading

# Specify the templates folder explicitly
app = Flask(__name__, template_folder="../templates")

# Global variables to store console output and alignment results
console_output = []
alignment_results = None

def perform_alignment(reference_file_path, query_file_path):
    """Perform alignment and store results in global variables."""
    global console_output, alignment_results
    console_output = []  # Reset console output

    # Load sequences
    console_output.append("Loading sequences...")
    reference_sequences = load_sequences(reference_file_path)
    query_sequences = load_sequences(query_file_path)

    # Perform alignment
    console_output.append("Performing alignment...")
    def log_to_console(message):
        console_output.append(message)

    best_match, best_score, probabilities, scores = find_closest_match(
        reference_sequences, query_sequences, log_callback=log_to_console
    )

    # Prepare results
    sorted_probabilities = sorted(probabilities.items(), key=lambda x: x[1], reverse=True)
    results = [
        {"breed": breed, "probability": f"{prob:.4f}", "score": scores[breed]}
        for breed, prob in sorted_probabilities[:10]  # Include top 10 entries
    ]

    console_output.append("Alignment complete.")
    alignment_results = (best_match, best_score, results)

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        # Get file paths from the form
        reference_file_path = request.form.get("reference_file", "data/dog_breeds.fa")
        query_file_path = request.form.get("query_file", "data/mystery.fa")

        # Redirect to loading page
        return redirect(url_for("loading", reference_file=reference_file_path, query_file=query_file_path))

    return render_template("index.html")

@app.route("/loading")
def loading():
    # Get file paths from query parameters
    reference_file_path = request.args.get("reference_file", "data/dog_breeds.fa")
    query_file_path = request.args.get("query_file", "data/mystery.fa")

    # Run alignment in a separate thread
    thread = threading.Thread(target=perform_alignment, args=(reference_file_path, query_file_path))
    thread.start()

    return render_template("loading.html")

@app.route("/console_output")
def console_output_route():
    """Serve console output as JSON."""
    global console_output
    complete = "Alignment complete." in console_output
    return jsonify({"output": console_output, "complete": complete})

@app.route("/results")
def results():
    # Retrieve alignment results
    global alignment_results
    best_match, best_score, results = alignment_results
    return render_template("results.html", results=results, best_match=best_match[0], best_score=best_score)

if __name__ == "__main__":
    app.run(debug=True)
