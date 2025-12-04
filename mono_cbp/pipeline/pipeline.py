"""Main pipeline orchestration for mono-cbp.

This module provides an interface to run the complete
CBP detection pipeline.
"""

import logging
from ..eclipse_masking import EclipseMasker
from ..transit_finding import TransitFinder
from ..model_comparison import ModelComparator
from ..injection_retrieval import TransitInjector
from ..utils import load_catalogue
from ..config import get_default_config, merge_config

logger = logging.getLogger('mono_cbp.pipeline')


class MonoCBPPipeline:
    """Main pipeline for circumbinary planet detection.

    This class orchestrates the complete workflow:
    1. Eclipse masking
    2. Transit finding
    3. Model comparison (vetting)
    4. (Optional) Injection-retrieval testing

    Attributes:
        catalogue: EB catalogue DataFrame
        config (dict): Pipeline configuration
        data_dir (str): Data directory path
        output_dir (str): Output directory path
        eclipse_masker (EclipseMasker): Eclipse masking component
        transit_finder (TransitFinder): Transit finding component
        model_comparator (ModelComparator): Model comparison component
        transit_injector (TransitInjector): Injection-retrieval component
        results (dict): Dictionary to store results from each pipeline step
    """

    def __init__(self, catalogue_path, data_dir='./data', output_dir='./results',
                 sector_times_path='../../catalogue/sector_times.csv', TEBC=False, transit_models_path=None, config=None):
        """Initialize the MonoCBPPipeline.

        Args:
            catalogue_path (str): Path to catalogue CSV file containing EB properties.
                Must contain columns: tess_id, period, bjd0, sectors, prim_pos, prim_width,
                sec_pos, sec_width.
            data_dir (str, optional): Directory containing light curve files. Defaults to './data'.
            output_dir (str, optional): Directory for output files. Defaults to './results'.
            sector_times_path (str, optional): Path to sector times CSV for Skye metric
            transit_models_path (str, optional): Path to transit models for injection-retrieval
            config (dict, optional): Configuration dictionary. Uses defaults if None.
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.config = merge_config(config, get_default_config()) if config else get_default_config()

        # Create output directory if it doesn't exist
        import os
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logger.info(f"Created output directory: {output_dir}")

        # Load catalogue
        self.catalogue = load_catalogue(catalogue_path, TEBC=TEBC)
        logger.info(f"Initialized pipeline with {len(self.catalogue)} targets")

        # Initialize components
        self.eclipse_masker = EclipseMasker(self.catalogue, data_dir=data_dir)

        self.transit_finder = TransitFinder(
            catalogue=self.catalogue,
            sector_times=sector_times_path,
            config=self.config
        )

        self.model_comparator = ModelComparator(config=self.config)

        # Injection-retrieval (optional)
        if transit_models_path:
            self.transit_injector = TransitInjector(
                transit_models_path=transit_models_path,
                catalogue=self.catalogue,
                config=self.config
            )
        else:
            self.transit_injector = None
            logger.info("Transit models not provided - injection-retrieval will not be available")

        self.results = {}

    def run(self, find_transits=True, vet_candidates=True,
            injection_retrieval=False, **kwargs):
        """Run the complete pipeline.

        Eclipse masking is always performed before transit finding to ensure
        eclipses are properly excluded from the transit search.

        Args:
            find_transits (bool, optional): Run transit finding. Defaults to True.
            vet_candidates (bool, optional): Run model comparison vetting. Defaults to True.
            injection_retrieval (bool, optional): Run injection-retrieval. Defaults to False.
            **kwargs: Additional keyword arguments passed to individual steps:
                - mask_eclipses_kwargs: Arguments for eclipse masking (e.g., force=True)
                - find_transits_kwargs: Arguments for transit finding
                - vet_candidates_kwargs: Arguments for vetting
                - injection_retrieval_kwargs: Arguments for injection-retrieval

        Returns:
            dict: Dictionary containing results from each pipeline step

        Example:
            >>> results = pipeline.run(
            ...     find_transits=True,
            ...     vet_candidates=True,
            ...     mask_eclipses_kwargs={'force': True},
            ...     find_transits_kwargs={'plot_output_dir': './plots'}
            ... )
        """
        logger.info("Starting mono-cbp pipeline")

        # Step 1: Eclipse masking (always run)
        logger.info("Step 1: Eclipse masking")
        self.mask_eclipses(**kwargs.get('mask_eclipses_kwargs', {}))

        # Step 2: Transit finding
        if find_transits:
            logger.info("Step 2: Transit finding")
            transit_results = self.find_transits(**kwargs.get('find_transits_kwargs', {}))
            self.results['transit_finding'] = transit_results

        # Step 3: Model comparison vetting
        if vet_candidates and find_transits:
            logger.info("Step 3: Model comparison vetting")
            # Get event snippets from transit finder if available
            event_snippets = self.transit_finder.results.get('event_snippets', [])

            # Prepare vetting kwargs
            vet_kwargs = kwargs.get('vet_candidates_kwargs', {})

            if len(event_snippets) > 0 and 'event_snippets' not in vet_kwargs:
                # Pass snippets directly to vetting
                logger.info(f"Passing {len(event_snippets)} event snippets")
                vet_kwargs['event_snippets'] = event_snippets
            elif len(event_snippets) == 0:
                # Fall back to file-based if no snippets in memory
                logger.info("No event snippets in memory, using file-based processing")

            vetting_results = self.vet_candidates(**vet_kwargs)
            self.results['vetting'] = vetting_results

        # Step 4 (optional): Injection-retrieval
        if injection_retrieval:
            if self.transit_injector is None:
                logger.warning("Injection-retrieval requested but transit models not provided - skipping")
            else:
                logger.info("Step 4: Injection-retrieval testing")
                inj_ret_results = self.run_injection_retrieval(**kwargs.get('injection_retrieval_kwargs', {}))
                self.results['injection_retrieval'] = inj_ret_results

        logger.info("Pipeline complete")
        return self.results

    def mask_eclipses(self, **kwargs):
        """Run eclipse masking on all files in data directory.

        Args:
            **kwargs: Passed to EclipseMasker.mask_all()
        """
        self.eclipse_masker.mask_all(**kwargs)
        logger.info("Eclipse masking complete")

    def find_transits(self, output_file='transit_events.txt', output_dir=None,
                     plot_output_dir=None):
        """Run transit finding on all files in data directory.

        Args:
            output_file (str, optional): Output filename. Defaults to 'transit_events.txt'.
            output_dir (str, optional): Output directory. Defaults to pipeline's output_dir.
            plot_output_dir (str, optional): Plot output directory. Defaults to None.

        Returns:
            pd.DataFrame: Detected transit events
        """
        if output_dir is None:
            output_dir = self.output_dir

        results = self.transit_finder.process_directory(
            self.data_dir,
            output_file=output_file,
            output_dir=output_dir,
            plot_output_dir=plot_output_dir
        )
        logger.info(f"Transit finding complete: {len(results)} events detected")
        return results

    def vet_candidates(self, event_snippets=None, event_snippets_dir=None,
                      output_file='vetting_results.csv', output_dir=None):
        """Run model comparison vetting on detected events.

        Args:
            event_snippets (list, optional): List of event snippet dictionaries to process in-memory.
                                            If provided, takes precedence over event_snippets_dir.
            event_snippets_dir (str, optional): Directory with event snippet .npz files.
                                               Defaults to data_dir/event_snippets if event_snippets is None.
            output_file (str, optional): Output filename. Defaults to 'vetting_results.csv'.
            output_dir (str, optional): Output directory. Defaults to pipeline's output_dir.

        Returns:
            pd.DataFrame: Vetting results with classifications
        """
        if output_dir is None:
            output_dir = self.output_dir

        # Use in-memory snippets if provided
        if event_snippets is not None:
            logger.info(f"Vetting {len(event_snippets)} candidates")
            results = self.model_comparator.compare_events(
                event_snippets,
                output_file=output_file,
                output_dir=output_dir
            )
        else:
            # Fall back to file-based processing
            if event_snippets_dir is None:
                import os
                event_snippets_dir = os.path.join(self.output_dir, 'event_snippets')

            logger.info(f"Vetting candidates from {event_snippets_dir}")
            results = self.model_comparator.compare_events(
                event_snippets_dir,
                output_file=output_file,
                output_dir=output_dir
            )

        logger.info(f"Model comparison complete: {len(results)} events vetted")
        return results

    def run_injection_retrieval(self, n_injections=None, output_file='injection_results.csv',
                               output_dir=None):
        """Run injection-retrieval testing for all transit models.

        Tests each transit model in transit_models.npz by injecting it into n_injections
        randomly selected light curves. The total number of tests will be
        n_injections * number_of_models.

        Args:
            n_injections (int, optional): Number of injection-retrieval tests to perform
                per transit model. If there are fewer files than requested injections,
                files will be randomly sampled with replacement. Defaults to config value.
            output_file (str, optional): Output filename. Defaults to 'injection_results.csv'.
            output_dir (str, optional): Output directory. Defaults to data_dir.

        Returns:
            pd.DataFrame: Injection-retrieval results with one row per injection test
        """
        if self.transit_injector is None:
            raise ValueError("Transit injector not initialized - provide transit_models_path")
        
        n_injections = self.config['injection_retrieval']['n_injections'] if n_injections is None else n_injections

        results = self.transit_injector.run_injection_retrieval(
            self.data_dir,
            n_injections=n_injections,
            output_file=output_file,
            output_dir=output_dir or self.data_dir
        )
        logger.info(f"Injection-retrieval complete")
        return results

    def plot_bin_phase_fold(self, tic_id, save_fig=False, save_path='.'):
        """Plot phase-folded light curve with eclipse masks.

        Args:
            tic_id (int): TIC ID to plot
            save_fig (bool, optional): Save figure. Defaults to False.
            save_path (str, optional): Save directory. Defaults to '.'.

        Returns:
            matplotlib.figure.Figure: Figure object
        """
        return self.eclipse_masker.plot_bin_phase_fold(tic_id, save_fig=save_fig, save_path=save_path)

    def plot_events(self, tic_id, event_number=None, save_fig=False, save_path='.', figsize=(12, 4)):
        """Plot detected transit events for a given TIC ID.

        Shows the detrended light curve with detected events highlighted, along with
        event properties (depth, width, SNR, phase). Can plot all events or a specific event.

        Args:
            tic_id (int): TIC ID to plot events for
            event_number (int, optional): Specific event number to plot (1-indexed).
                                         If None, plots all events for this TIC. Defaults to None.
            save_fig (bool, optional): Whether to save the figure. Defaults to False.
            save_path (str, optional): Directory to save figure. Defaults to '.'.
            figsize (tuple, optional): Figure size (width, height). Defaults to (12, 4) per event.

        Returns:
            None

        Raises:
            ValueError: If no events found for the given TIC ID or event number
            RuntimeError: If transit finding has not been run yet

        Example:
            >>> # Plot all events for TIC 260128333
            >>> pipeline.plot_events(260128333)
            >>>
            >>> # Plot only the first event
            >>> pipeline.plot_events(260128333, event_number=1)
        """
        import matplotlib.pyplot as plt
        import numpy as np
        import os

        # Check if transit finding has been run
        if not hasattr(self.transit_finder, 'results') or len(self.transit_finder.results.get('tics', [])) == 0:
            raise RuntimeError("No transit events found. Run pipeline.find_transits() or pipeline.run() first.")

        # Get all events for this TIC
        tic_str = str(tic_id)
        event_indices = [i for i, tic in enumerate(self.transit_finder.results['tics']) if tic == tic_str]

        if len(event_indices) == 0:
            raise ValueError(f"No events found for TIC {tic_id}")

        # Filter to specific event if requested
        if event_number is not None:
            if event_number < 1 or event_number > len(event_indices):
                raise ValueError(f"Event number {event_number} out of range. TIC {tic_id} has {len(event_indices)} event(s).")
            event_indices = [event_indices[event_number - 1]]

        # Get event snippets if available
        event_snippets = self.transit_finder.results.get('event_snippets', [])

        # Find matching snippets for this TIC
        matching_snippets = []
        for idx in event_indices:
            # Try to find the corresponding snippet
            if len(event_snippets) > idx:
                snippet = event_snippets[idx]
                if snippet['tic'] == tic_id:
                    matching_snippets.append((idx, snippet))

        if len(matching_snippets) == 0:
            logger.warning(f"No event snippet data available for TIC {tic_id}. "
                          "Event snippets may not have been generated during transit finding. "
                          "Set 'generate_event_snippets': True in config.")
            return

        # Create figure with subplots for each event
        n_events = len(matching_snippets)
        fig, axes = plt.subplots(n_events, 1, figsize=(figsize[0], figsize[1] * n_events),
                                squeeze=False, constrained_layout=True)

        for plot_idx, (event_idx, snippet) in enumerate(matching_snippets):
            ax = axes[plot_idx, 0]

            # Extract data from snippet
            time = snippet['time']
            flux = snippet['flux']
            flux_err = snippet['flux_err']
            event_time = snippet['event_time']
            event_width = snippet['event_width']

            # Get event properties from results
            sector = self.transit_finder.results['sectors'][event_idx]
            depth = self.transit_finder.results['event_depths'][event_idx]
            duration = self.transit_finder.results['event_durations'][event_idx]
            snr = self.transit_finder.results['event_snrs'][event_idx]
            phase = self.transit_finder.results['event_phases'][event_idx]

            # Plot light curve with error bars (straight lines, no caps)
            ax.errorbar(time, flux, yerr=flux_err, fmt='o', color='navy',
                       alpha=0.6, markersize=4, capsize=0)

            # Highlight the event region
            event_start = event_time - event_width / 2
            event_end = event_time + event_width / 2
            ax.axvspan(event_start, event_end, alpha=0.2, color='red')
            ax.axvline(event_time, color='red', linestyle='--', alpha=0.7, linewidth=1.5)

            # Add horizontal line at flux=1
            ax.axhline(1.0, color='gray', linestyle=':', alpha=0.5, linewidth=1)

            # Format plot
            ax.set_xlabel('Time (BJD - 2457000)', fontsize=11)
            ax.set_ylabel('Normalized Flux', fontsize=11)

            # Create title with event info
            event_num = event_idx - event_indices[0] + 1 if event_number is None else event_number
            title = (f'TIC {tic_id} - Event {event_num} (Sector {sector})\n'
                    f'Depth: {depth:.4f} | Duration: {duration:.3f} d | SNR: {snr:.1f} | Phase: {phase:.3f}')
            ax.set_title(title, fontsize=12, fontweight='bold')

            ax.grid(True, alpha=0.3, linestyle='--')

            # Set y-axis limits with some padding
            flux_range = np.max(flux) - np.min(flux)
            ax.set_ylim(np.min(flux) - 0.1 * flux_range, np.max(flux) + 0.1 * flux_range)

        if save_fig:
            if event_number is not None:
                filename = f"TIC_{tic_id}_event_{event_number}.png"
            else:
                filename = f"TIC_{tic_id}_all_events.png"
            output_path = os.path.join(save_path, filename)
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Saved event plot to {output_path}")
            plt.close(fig)
        else:
            plt.show()
