import RNA
import numpy as np

def find_triggers(sequence, min_length=10, stride=5):
    """
    Find potential trigger sequences in a given RNA sequence.

    Parameters:
    sequence (str): The RNA sequence to analyze.
    min_length (int): The minimum length of the trigger sequences to consider.
    stride (int): The step size for the sliding window approach.

    Returns:
    list: A list of potential trigger sequences.
    """

    #Alright we're just gonna go with a sliding window approach
    
    #calculate number of windows (last window may be smaller than min_length)
    num_windows = (len(sequence) - min_length) // stride
    MFEs = []
    for i in range(num_windows + 1):
        start = i * stride
        end = start + min_length
        window_seq = sequence[start:end]
        
        # Calculate the minimum free energy (MFE) of the window sequence
        mfe = RNA.fold(window_seq)[0]  # Get the MFE value
        MFEs.append((window_seq, mfe))

    #Add a flat bonus to the MFE values near the edges of the sequence to encourage triggers near the ends
    #something like a very thick Gaussian centered at the edges of the sequence
    for i in range(len(MFEs)):
        gaussian_bonus = np.exp(-(len(sequence)/2 - i * stride) ** 2) / (2 * (len(sequence) / 10) ** 2))
        MFEs[i] = (MFEs[i][0], MFEs[i][1] - gaussian_bonus)

    # Sort the sequences by their MFE values (higher is better)
    MFEs.sort(key=lambda x: x[1], reverse=True)
    return MFEs


    

        