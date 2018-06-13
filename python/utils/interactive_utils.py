#!/usr/bin/env python
"""
Module containing some classes and functions that help with interactive things
in jupyter notebooks
"""

import ipywidgets as widgets

from collections import OrderedDict

from utils.misc_helpers import log_key_error


class InteractiveSelection(object):
    """
    Class helping with choices by holding all the necessary state
    """
    def __init__(self, choices=None, select=False):
        """
        Create, possibly with a predefined set of choices
        If select is True all possibilities will be selected in the beginning
        """
        self.choices = OrderedDict()
        self.curr_selected = OrderedDict()
        self.sel_boxes = []
        if choices is not None:
            for key, val in choices.iteritems():
                self.choices[key] = val
                self.curr_selected[key] = select


    def add_choice(self, name, value, select=False):
        """Add another choice with a name and the value"""
        self.choices[name] = value
        self.curr_selected[name] = select


    @log_key_error
    def rm_choice(self, name):
        """Remove a choice with a given name"""
        del self.choices[name]
        del self.curr_selected[name]


    def get_selected(self):
        """Get all the currently selected choices"""
        sel_choices = OrderedDict()
        for key, selected in self.curr_selected.iteritems():
            if selected:
                sel_choices[key] = self.choices[key]
        return sel_choices


    def get_selected_values(self):
        """Get a list of the selected values"""
        return [v for (k, v) in self.choices.iteritems()
                if self.curr_selected[k]]

    def get_selected_keys(self):
        """Get a list of the selected keys"""
        return [k for k in self.choices if self.curr_selected[k]]


    def select(self, name):
        """Select the choice with the given name"""
        self._set_selected(name, True)


    def deselect(self, name):
        """Deselect the choice with the given name"""
        self._set_selected(name, False)


    def select_all(self):
        """Select all elements"""
        for choice in self.curr_selected:
            self.curr_selected[choice] = True


    def deselect_all(self):
        """Deselect all elements"""
        for choice in self.curr_selected:
            self.curr_selected[choice] = False


    def interactive_selection(self):
        """Get an interactive selection"""
        for choice in self.choices:
            self.sel_boxes.append(widgets.Checkbox(self.curr_selected[choice],
                                                   description=choice))
        widgets.interact_manual(self.update_selection,
                                **{s.description: s.value for s in
                                   self.sel_boxes})
        widgets.interact_manual(self.get_selected)
        widgets.interact_manual(self.select_all)
        widgets.interact_manual(self.deselect_all)


    def update_selection(self, **kwargs):
        """bulk update to a number of selections"""
        for choice, selected in kwargs.iteritems():
            self._set_selected(choice, selected)


    @log_key_error
    def _set_selected(self, name, selected):
        """Helper function for setting the status of a given choice"""
        # Do not set non-existing values
        if not name in self.curr_selected:
            raise KeyError
        self.curr_selected[name] = selected
