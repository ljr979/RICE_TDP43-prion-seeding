
"""Generate Cellpose masks from microscopy TIFF images.

The script discovers TIFF files in a configured input directory, adjusts image
intensity, runs the selected Cellpose model, and saves mask TIFFs and optional
quality-control plots.
"""

# =====================
# CONFIGURATION SECTION
# =====================
#important - toggle this to FALSE if you want to run this on a big dataset/a full dataset. This will STOP all of the plotting of every image (this slows down your computer) I would leave it as TRUE when doing test runs to check parameters at the start of each experiment.
plot_all_figures = False 
#-----------------
#start here by looping over all files in one treatment, then we can extend to more treatments later
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from cellpose import models
from cellpose.plot import mask_overlay
import matplotlib
if not plot_all_figures:
    matplotlib.use('Agg') 
from matplotlib import colors
import random
from skimage import exposure

def plot_cellpose_results(img, masks):
    """Display an image, its instance labels, and their overlay.

    Args:
        img: Raw two-dimensional image or multichannel image array.
        masks: Cellpose label array with the same spatial dimensions.

    Returns:
        matplotlib.figure.Figure: The three-panel quality-control figure.
    """

    # Select the first channel for display when the input is multichannel.
    if img.ndim == 3 and img.shape[-1] > 1:
        disp = img[..., 0]
    else:
        disp = img

    # Scale image intensities for display.
    disp_norm = exposure.rescale_intensity(disp, out_range=(0, 1))

    # Assign a random color to each mask label.
    cmap = colors.ListedColormap(np.random.rand(256, 3))

    # Create an overlay of the image and segmentation labels.
    overlay = mask_overlay(disp_norm, masks)

    # Arrange the source image, labels, and overlay in one figure.
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].imshow(disp_norm, cmap="gray")
    axes[0].set_title("Raw image")
    axes[0].axis("off")

    axes[1].imshow(masks, cmap=cmap, interpolation="nearest")
    axes[1].set_title("Masks (labels)")
    axes[1].axis("off")

    axes[2].imshow(overlay)
    axes[2].set_title("Overlay")
    axes[2].axis("off")

    plt.tight_layout()

    return fig

def run_cellpose(extracted_part, model_used='cyto2'):
    """Run the selected Cellpose model and return its segmentation labels.

    Args:
        extracted_part: Image array supplied to Cellpose.
        model_used: Cellpose model name.

    Returns:
        numpy.ndarray: Segmentation labels produced by Cellpose.
    """
    model = models.Cellpose(gpu=False, model_type=model_used)
    masks, _, _, _ = model.eval(
                extracted_part,
                diameter=0,
                flow_threshold=0.9,
                cellprob_threshold=-5,
                anisotropy=1.0,
                do_3D=False,
                stitch_threshold=0.0,
                channels=[0, 1],
            )
    
    return masks

def adjust_image_for_cellpose(plot_all_figures, normalize_image, contrast_stretch, gamma_adjust, filename, image):
    """Normalize, contrast-adjust, and convert an image for Cellpose input.

    Args:
        plot_all_figures: Whether to display the adjusted image.
        normalize_image: Image normalization function.
        contrast_stretch: Percentile-based contrast adjustment function.
        gamma_adjust: Gamma correction function.
        filename: Source image name used in the display title.
        image: Raw image array.

    Returns:
        numpy.ndarray: Adjusted 8-bit image array.
    """
    # Prepare a single adjusted image for segmentation.
    gfp_norm = normalize_image(image)


    # Apply mild contrast stretching between the selected percentiles.
    gfp_adj = contrast_stretch(gfp_norm, low=1, high=99)


    # Brighten midtones using the configured gamma value.
    gfp_adj = gamma_adjust(gfp_adj, gamma=0.8)


    # Use the adjusted image as the Cellpose input.
    combined = gfp_adj 

    if plot_all_figures:
        rgb = np.zeros((gfp_adj.shape[0], gfp_adj.shape[1], 3)) 
        rgb[..., 1] = gfp_adj
        plt.figure(figsize=(6, 6))
        plt.imshow(rgb, origin='lower')
        plt.title("False-colour overlay (red=ch1, green=ch0)- weighted sum of channels")
        plt.axis('off')
        plt.show()
                

    combined_8bit = (combined / combined.max() * 255).astype(np.uint8)

    plt.figure(figsize=(6, 6))
    plt.imshow(combined_8bit, cmap='gray', origin='lower')
    plt.title(f'{filename}_adjusted')

    return combined_8bit

def normalize_image(img):
    """Scale image intensities to the interval from zero to one."""
    # Convert to floating-point values before normalization.
    img = img.astype(float)
    return (img - img.min()) / (img.max() - img.min() + 1e-8)

def contrast_stretch(img, low=1, high=99):
    """Stretch image contrast using lower and upper intensity percentiles."""
    p_low, p_high = np.percentile(img, (low, high))
    stretched = np.clip((img - p_low) / (p_high - p_low + 1e-8), 0, 1)
    return stretched

def gamma_adjust(img, gamma=1.0):
    """Apply gamma correction; values below one brighten midtones."""
    return np.power(img, gamma)

    

if __name__ == '__main__':

    # Set the microscopy image directory and optional alternate output folder.
    input_folder = 'Z:/Raw Data/Lauren/Experiment115_SOD1vars-ZMC1/rep_1/TileScan 1/'
    results_path = 'results/'
    # Select TIFF filenames containing this channel identifier.
    file_identifier = 'ch00'

    # Find TIFF images recursively beneath the configured input directory.
    filepaths = [[f'{root}/{name}' for name in files if '.tif' in name] for root, dirs, files in os.walk(f'{input_folder}/')]
    # Exclude masks already generated by a previous run.
    filepaths = [item for sublist in filepaths for item in sublist if f'{file_identifier}' in item and 'mask' not in item]

    # Optionally sample a small subset for a test run.
    test = False
    # Set the number of images in the test subset.
    num_tests = 1
    # Control whether test-run overlays are saved.
    save_figs = False
    # Select the Cellpose model.
    model_used = 'cyto2'
    # Choose whether masks are written beside the source images or in results_path.
    diff_output_folder = False

    if test:
        filepaths = random.sample(filepaths, num_tests)



    # Process each selected image.
    for image_path in filepaths:
        # Normalize path separators before parsing the filename.
        image_path = image_path.replace('\\', '/')
        image_path = image_path.replace('//', '/')
        print(f'currently analysing {image_path}')
        # Derive the source image stem for its output filenames.
        filename = image_path.split('/')[-1].split('.')[0]

        # Write outputs beside the source image by default.
        save_path = f'{"/".join(image_path.split("/")[0:(len(image_path.split("/"))-1)])}/'

        if diff_output_folder:
            save_path = f'{results_path}/'
        # Load the source TIFF image.
        image = tifffile.imread(image_path)

        if test:
            # Display the original grayscale image during a test run.
            fig, axes = plt.subplots()
            plt.imshow(image, cmap="gray")
            plt.title(f'{filename}')

        # Adjust image intensity before segmentation.
        combined_8bit = adjust_image_for_cellpose(plot_all_figures, normalize_image, contrast_stretch, gamma_adjust, filename, image)
        # Segment the adjusted image with Cellpose.
        masks = run_cellpose(combined_8bit, model_used)

        # Optionally inspect and save a quality-control overlay.
        if test:
            plot_cellpose_results(combined_8bit, masks)
            if save_figs:
                plt.savefig(f'{save_path}{filename}_mask_overlay.png')
            plt.show()

        tifffile.imwrite(f'{save_path}{filename}_mask.tif', masks.astype(np.uint16))
        
        plt.close()